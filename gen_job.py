import argparse
import subprocess
import re
from pathlib import Path

_RE_HW = re.compile(r'hw(\d+)')
_RE_E2MAX = re.compile(r'_emax(\d+)_e2max(\d+)\.')
_RE_E3MAX = re.compile(r'_emax(\d+)_e2max(\d+)_e3max(\d+)\.')

REPO_ROOT = Path(__file__).resolve().parent

partition = "c128m1024"
cpus = 128

exe = str(REPO_ROOT / "src" / "imsrg++")
params = {}

flag = "EM1.8_2.0"
params["fmt2"] = "me2j"
dir_nn = "/tns/public/Forces/EM1.8_2.0/2BME/"
params["2bme"] = dir_nn + "TwBME_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz"
params["fmt3"] = ""
params["3bme_type"] = "no2b"
params["no2b_precision"] = "single"
dir_3n = "/tns/public/Forces/EM1.8_2.0/3BME/"
params["3bme"] = dir_3n + "ThBME_NO2B_EM1.8_2.0_hw16_emax14_e2max28_e3max18.me3j.gz"
params["emax"] = 4
params["e3max"] = 4

# For a packed normal-ordered Hamiltonian produced by the normal-order program:
#   params["fmt2"] = "no2bpack"
#   params["2bme"] = "<normal-order output .bin>"
#   params["3bme"] = "none"
# The packed file already contains the normal-ordered 0B/1B/2B pieces, including
# the 3N contribution from the generation step. Keep params["emax"] and
# params["e3max"] consistent with the truncations used to generate the pack.


params["reference"] = "O18"
params["valence_space"] = "sd-shell"
params["Operators"] = ""

core = params["valence_space"].split(",")[0]

params["basis"] = "HF"
params["method"] = "magnus"


params["core_generator"] = "atan"
params["valence_generator"] = "shell-model-atan"


def check_and_make_dir(dir: str):
    path = Path(dir)
    if not path.exists():
        path.mkdir(parents=True, exist_ok=True)


def extract_emax_e2max(file: str):
    match = re.search(_RE_E2MAX, file)
    if match:
        return tuple(map(int, match.groups()))
    else:
        print("not match emax e2max")
        exit(-1)


def extract_emax_e2max_e3max(file3: str):
    if file3 == "none":
        return (-1, -1, -1)
    match = re.search(_RE_E3MAX, file3)
    if match:
        return tuple(map(int, match.groups()))
    else:
        print("not match emax e2max e3max")
        exit(-1)


def extract_hw(file, file3):
    match1 = re.search(_RE_HW, file)
    hw1 = int(match1.group(1))
    if file3 == "none":
        return hw1
    match2 = re.search(_RE_HW, file3)
    if match2:
        hw2 = int(match2.group(1))

        if hw1 == hw2:
            return hw1
        else:
            print("hw of two body interaction and three body interaction are different")
            exit(-1)


def add_header(output, partition=partition, cpus=cpus):
    lib_dir = REPO_ROOT / "src"
    header_list = [
        "#!/bin/bash -l",
        f"#SBATCH --partition={partition}",
        "#SBATCH -J IMSRG",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks-per-node=1",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH -o {output}",
        "set -e",
        f'cd "{REPO_ROOT}"',
        "source ./sourceme.sh",
        f"export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-{cpus}}}",
        f'export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:{lib_dir}"',
    ]
    header = '\n'.join(header_list)
    return header


def generate_slurm(params: dict):
    file = params["2bme"]
    file3 = params["3bme"]
    hw = extract_hw(file, file3)
    params["hw"] = hw
    A = int(re.findall(r'(\d+)', params["reference"])[0])
    params["A"] = A

    emax_nn, e2max_nn = extract_emax_e2max(file)
    emax_3n, e2max3_3n, e3max_3n = extract_emax_e2max_e3max(file3)

    params["file2e1max"] = emax_nn
    params["file2e2max"] = e2max_nn

    params["file3e1max"] = emax_3n
    params["file3e2max"] = e2max3_3n
    params["file3e3max"] = e3max_3n

    reference = params["reference"]
    emax = params["emax"]
    e3max = params["e3max"]
    prefix = f"{flag}_{core.lower()}_{reference.lower()}_hw{hw}_emax{emax}_e3max{e3max}"

    result_dir = REPO_ROOT / "result" / flag / f"{core.lower()}_{reference.lower()}_hw{hw}_emax{emax}_e3max{e3max}"
    script_file = result_dir / f"run_{prefix}.sh"
    scratch_dir = result_dir / "scratch"

    check_and_make_dir(result_dir)
    check_and_make_dir(scratch_dir)

    params["flowfile"] = str(result_dir / f"{prefix}.dat")
    params["intfile"] = str(result_dir / f"{prefix}")

    params_str = " ".join(f"{key}={val}" for key, val in params.items())

    with open(script_file, 'w') as f:
        f.write(add_header(result_dir / f"log_{prefix}_%j_%N.txt"))
        f.write('\n')
        f.write(f"{exe} {params_str}\n")
    return script_file


def run_smoke_test():
    result_dir = REPO_ROOT / "result" / "smoke_test"
    check_and_make_dir(result_dir)
    script_file = result_dir / "run_help.sh"
    log_file = result_dir / "help.log"

    with open(script_file, "w") as f:
        f.write(add_header(log_file))
        f.write("\n")
        f.write(f'{exe} help > "{log_file}" 2>&1\n')

    subprocess.run(["bash", str(script_file)], check=True)
    print(f"smoke test passed: {log_file}")


def parse_args():
    parser = argparse.ArgumentParser(description="Generate and optionally run an IMSRG Slurm job script.")
    parser.add_argument("--generate-only", action="store_true", help="only generate the Slurm script")
    parser.add_argument("--submit", action="store_true", help="submit the generated script with sbatch")
    parser.add_argument("--smoke-test", action="store_true", help="run a lightweight script that calls imsrg++ help")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    check_and_make_dir(REPO_ROOT / "result")

    if args.smoke_test:
        run_smoke_test()
        exit(0)

    script_file = generate_slurm(params)

    if args.generate_only:
        print(script_file)
        exit(0)

    if args.submit:
        subprocess.run(
            ["sbatch", str(script_file)],
            capture_output=True,
            text=True,
            check=True,
        )
        print(f"sbatch {script_file}")
        exit(0)

    run_mode = input("submit job: y/n\n")
    if run_mode.lower()[0] == "y":
        subprocess.run(
            ["sbatch", str(script_file)],
            capture_output=True,
            text=True,
            check=True,
        )
        print(f"sbatch {script_file}")
    else:
        subprocess.run(
            ["bash", str(script_file)],
            check=True,
        )
