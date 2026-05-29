import argparse
import subprocess
import re
from pathlib import Path

_RE_HW = re.compile(r'hw(\d+)')
_RE_E2MAX = re.compile(r'_emax(\d+)_e2max(\d+)\.')
_RE_E3MAX = re.compile(r'_emax(\d+)_e2max(\d+)_e3max(\d+)\.')
_RE_NO2BPACK = re.compile(r'_hw(\d+)_emax(\d+)_e3max(\d+)\.')

REPO_ROOT = Path(__file__).resolve().parent

# point7 Slurm partitions:
#   c128m1024: 128 CPU cores, about 1 TB memory per node (default)
#   c128m512 : 128 CPU cores, about 512 GB memory per node
#   compute_C: 96 CPU cores per node
#   compute_A: 28 CPU cores per node
# point7 has no accounting/QOS configured (accounting_storage/none), so any
# --qos is silently ignored there; keep qos=low so the wm2 setting is unchanged.
partition = "c128m512"
qos = "low"
cpus = 64
# Optionally pin to a specific node, e.g. "node2", to keep this job off a busy
# node and on its own node (point7 does not track memory, so co-located jobs
# share RAM with no protection). Set to None / "" to let Slurm choose.
nodelist = "node3"

exe = str(REPO_ROOT / "build" / "imsrg++")
params = {}

flag = "EM1.8_2.0"
params["fmt2"] = "me2j"
params["2bme"] = "/tns/public/Forces/EM1.8_2.0/2BME/TwBME_N3LO_EM500_srg1.8_hw12_emax14_e2max28.me2j.gz"
params["fmt3"] = ""
params["3bme_type"] = "no2b"
params["no2b_precision"] = "single"
params["3bme"] = "/tns/public/Forces/EM1.8_2.0/3BME/ThBME_NO2B_EM1.8_2.0_hw12_emax14_e2max28_e3max24.me3j.gz"
params["emax"] = 14
params["e3max"] = 24
params["BetaCM"] = 0.0
params["denominator_delta"] = 0.0
params["denominator_delta_orbit"] = "all"
params["nucleon_mass_correction"] = "true"

# For a packed normal-ordered Hamiltonian produced by the normal-order program:
#   params["fmt2"] = "no2bpack"
#   params["2bme"] = "<normal-order output .bin>"
#   params["3bme"] = "none"
# The packed file already contains the normal-ordered 0B/1B/2B pieces, including
# the 3N contribution from the generation step. Keep params["emax"] and
# params["e3max"] consistent with the truncations used to generate the pack.


params["reference"] = "Fe57"

# Model-space naming convention:
# - params["valence_space"] is the short name used in paths and output prefixes.
#   If it is a built-in IMSRG space such as "p-shell", it is also parsed by imsrg++.
# - params["custom_valence_space"], when set, is the physical model-space
#   definition parsed by imsrg++. Its format is "<core>,<orbit>,...", for example
#   "He4,p0p3,n0p3,p0p1,n0p1,p0d5,n0d5,p1s1,n1s1".
# - Do not hand-write a short name and custom definition that describe different
#   spaces; verify the intended orbit content before generating/submitting.
params["valence_space"] = "fp-shell"
# fp-shell is a built-in IMSRG model space (Ca40 core + f7/2,p3/2,f5/2,p1/2 for
# protons and neutrons), so imsrg++ parses the short name directly. For a
# built-in space, comment out custom_valence_space instead of removing it.
# params["custom_valence_space"] = "Ca40,p0f7,n0f7,p1p3,n1p3,p0f5,n0f5,p1p1,n1p1"
params["Operators"] = "Sigma,SigmaTau3,Ltau3"

space_tag = params["valence_space"]

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


def extract_no2bpack_hw_emax_e3max(file: str):
    match = re.search(_RE_NO2BPACK, file)
    if match:
        return tuple(map(int, match.groups()))
    else:
        print("not match no2bpack hw emax e3max")
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


def path_token(value):
    return f"{value:g}"


def add_header(output, partition=partition, qos=qos, cpus=cpus, nodelist=nodelist):
    lib_dir = REPO_ROOT / "build"
    header_list = [
        "#!/bin/bash -l",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --qos={qos}",
        "#SBATCH -J IMSRG",
    ]
    if nodelist:
        header_list.append(f"#SBATCH --nodelist={nodelist}")
    header_list += [
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
    A = int(re.findall(r'(\d+)', params["reference"])[0])
    params["A"] = A

    if params["fmt2"] == "no2bpack":
        hw, emax_nn, e3max_nn = extract_no2bpack_hw_emax_e3max(file)
    else:
        hw = extract_hw(file, file3)
        emax_nn, e2max_nn = extract_emax_e2max(file)
        emax_3n, e2max3_3n, e3max_3n = extract_emax_e2max_e3max(file3)

        params["file2e1max"] = emax_nn
        params["file2e2max"] = e2max_nn
        params["file3e1max"] = emax_3n
        params["file3e2max"] = e2max3_3n
        params["file3e3max"] = e3max_3n

    params["hw"] = hw
    params["hwBetaCM"] = hw

    reference = params["reference"]
    emax = params["emax"]
    e3max = params["e3max"]
    lawson_tag = f"beta{path_token(params['BetaCM'])}"
    delta_tag = f"delta{path_token(params['denominator_delta'])}"
    run_tag = f"{lawson_tag}_{delta_tag}"
    prefix = f"{flag}_{space_tag.lower()}_{reference.lower()}_hw{hw}_emax{emax}_e3max{e3max}_{run_tag}"

    result_dir = REPO_ROOT / "result" / flag / f"{space_tag.lower()}_{reference.lower()}_hw{hw}_emax{emax}_e3max{e3max}_{run_tag}"
    script_file = result_dir / f"run_{prefix}.sh"
    scratch_dir = result_dir / "scratch"

    check_and_make_dir(result_dir)
    check_and_make_dir(scratch_dir)

    params["flowfile"] = str(result_dir / f"{prefix}.dat")
    params["intfile"] = str(result_dir / f"{prefix}")

    params_str = " ".join(f"{key}={val}" for key, val in params.items())

    with open(script_file, 'w') as f:
        f.write(add_header(result_dir / f"log_{prefix}_%j.txt"))
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
        result = subprocess.run(
            ["sbatch", str(script_file)],
            capture_output=True,
            text=True,
            check=True,
        )
        print(f"sbatch {script_file}")
        if result.stdout.strip():
            print(result.stdout.strip())
        exit(0)

    run_mode = input("submit job: y/n\n")
    if run_mode.lower()[0] == "y":
        result = subprocess.run(
            ["sbatch", str(script_file)],
            capture_output=True,
            text=True,
            check=True,
        )
        print(f"sbatch {script_file}")
        if result.stdout.strip():
            print(result.stdout.strip())
    else:
        subprocess.run(
            ["bash", str(script_file)],
            check=True,
        )
