import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
import subprocess
import re
import shlex
from pathlib import Path
from typing import Optional

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


@dataclass(frozen=True)
class MRJschemeSettings:
    nucleus: str
    nrefmax: int
    emax: int
    interaction: Path
    reference_file: Path
    target_s: float
    start_s: float = 0.0
    method: str = "flow"
    ds_0: float = 1e-2
    dsmax: float = 0.5
    ode_tolerance: float = 1e-8
    eta_criterion: float = 1e-6
    partition: str = "c128m512"
    qos: str = "low"
    cpus: int = 64
    nodelist: Optional[str] = None
    label: Optional[str] = None
    result_root: Path = REPO_ROOT / "result" / "mr-jscheme-flow"
    executable: Path = REPO_ROOT / "build" / "src" / "imsrg++"


@dataclass(frozen=True)
class MRJschemeInputSettings:
    nucleus: str
    nrefmax: int
    emax: int
    minipack: Path
    source_reference: Path
    partition: str = "c128m512"
    qos: str = "low"
    cpus: int = 64
    nodelist: Optional[str] = None
    label: Optional[str] = None
    result_root: Path = REPO_ROOT / "result" / "mr-jscheme-inputs"
    converter: Path = (
        REPO_ROOT / "prototype" / "mrimsrg" / "build" /
        "mrimsrg_minipack_to_j64"
    )
    pyimsrg_dir: Path = REPO_ROOT / "build" / "src"


@dataclass(frozen=True)
class MRNCSMReadbackSettings:
    no2bpack: Path
    proton_number: int
    neutron_number: int
    nmax: int
    states: int = 3
    max_iter: int = 300
    partition: str = "c128m512"
    qos: str = "low"
    cpus: int = 64
    nodelist: Optional[str] = None
    label: Optional[str] = None
    sample_rss: bool = True
    result_root: Path = REPO_ROOT / "result" / "mr-ncsm-readback"
    executable: Path = REPO_ROOT / "prototype" / "mrimsrg" / "build-emax6" / "mrimsrg_validate"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _repository_commit(repository: Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(repository), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    commit = result.stdout.strip()
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise RuntimeError("git returned an invalid repository commit")
    return commit


def _repository_state_files(repository: Path) -> tuple[tuple[Path, str], ...]:
    """Freeze files resolving HEAD without requiring Git on a compute node.

    point7 login nodes provide ``git``, while some compute-node images do not.
    Hashing the worktree HEAD file and, for an attached checkout, its loose ref
    (or packed-refs fallback) preserves the runtime checkout guard without
    depending on a compute-node Git installation.
    """

    def git_path(name: str) -> Path:
        result = subprocess.run(
            ["git", "-C", str(repository), "rev-parse", "--git-path", name],
            check=True,
            capture_output=True,
            text=True,
        )
        path = Path(result.stdout.strip())
        if not path.is_absolute():
            path = repository / path
        return path.resolve()

    head = git_path("HEAD")
    if not head.is_file():
        raise FileNotFoundError(f"repository HEAD file is unavailable: {head}")
    paths = [head]
    head_text = head.read_text(encoding="ascii").strip()
    if head_text.startswith("ref: "):
        ref_name = head_text.removeprefix("ref: ")
        if not re.fullmatch(r"refs/[A-Za-z0-9._/-]+", ref_name):
            raise RuntimeError(
                f"repository HEAD contains an invalid ref: {ref_name}"
            )
        ref_file = git_path(ref_name)
        if ref_file.is_file():
            paths.append(ref_file)
        else:
            packed_refs = git_path("packed-refs")
            if not packed_refs.is_file():
                raise FileNotFoundError(
                    f"repository ref is neither loose nor packed: {ref_name}"
                )
            paths.append(packed_refs)
    return tuple((path, _sha256(path)) for path in paths)


def _repository_state_checks(
    state_files: tuple[tuple[Path, str], ...],
) -> str:
    return "\n".join(
        f"echo {shlex.quote(digest + '  ' + str(path))} | sha256sum -c -"
        for path, digest in state_files
    )


def _validate_mr_jscheme_settings(settings: MRJschemeSettings) -> None:
    if not re.fullmatch(r"(?:He4|Be8|C12|O16)", settings.nucleus):
        raise ValueError("MR J-scheme nucleus must be He4, Be8, C12, or O16")
    allowed_nrefmax = {"He4": (0, 2), "Be8": (0,), "C12": (0,), "O16": (0, 2)}
    if settings.nrefmax not in allowed_nrefmax[settings.nucleus]:
        raise ValueError("unsupported nucleus/Nrefmax MR reference")
    if settings.emax < 2 or settings.emax % 2:
        raise ValueError("MR J-scheme emax must be an even integer >= 2")
    if settings.method not in ("flow", "flow_adaptive", "flow_RK4"):
        raise ValueError("MR J-scheme production currently supports direct flow only")
    if not settings.start_s >= 0.0 or not settings.target_s > settings.start_s:
        raise ValueError("MR J-scheme requires 0 <= start_s < target_s")
    if min(settings.ds_0, settings.dsmax, settings.ode_tolerance,
           settings.eta_criterion) <= 0.0:
        raise ValueError("MR flow steps, tolerance, and eta criterion must be positive")
    if settings.cpus < 1:
        raise ValueError("MR J-scheme cpus must be positive")
    if settings.partition not in ("c128m1024", "c128m512", "compute_C", "compute_A"):
        raise ValueError("unsupported point7 partition")
    if settings.qos != "low":
        raise ValueError("MR J-scheme production requires qos=low")
    for value, field in ((settings.nodelist, "nodelist"), (settings.label, "label")):
        if value is not None and not re.fullmatch(r"[A-Za-z0-9_-]+", value):
            raise ValueError(f"{field} contains unsupported characters")
    if not settings.interaction.is_file():
        raise FileNotFoundError(f"MR J-scheme interaction is unavailable: {settings.interaction}")
    if not settings.reference_file.is_file():
        raise FileNotFoundError(f"MR reference is unavailable: {settings.reference_file}")
    if not settings.executable.is_file():
        raise FileNotFoundError(f"imsrg++ executable is unavailable: {settings.executable}")
    if not (settings.executable.parent / "libIMSRG.so").is_file():
        raise FileNotFoundError("libIMSRG.so is unavailable beside imsrg++")
    if not (REPO_ROOT / "sourceme.sh").is_file():
        raise FileNotFoundError("repository sourceme.sh is unavailable")


def generate_mr_jscheme_slurm(settings: MRJschemeSettings) -> Path:
    """Generate one production C++ J-scheme MR direct-flow Slurm script."""
    _validate_mr_jscheme_settings(settings)
    interaction = settings.interaction.resolve()
    reference_file = settings.reference_file.resolve()
    executable_path = settings.executable.resolve()
    shared_library = (settings.executable.parent / "libIMSRG.so").resolve()
    environment_script = (REPO_ROOT / "sourceme.sh").resolve()
    result_root = settings.result_root.resolve()
    interaction_sha256 = _sha256(interaction)
    reference_sha256 = _sha256(reference_file)
    executable_sha256 = _sha256(executable_path)
    shared_library_sha256 = _sha256(shared_library)
    environment_script_sha256 = _sha256(environment_script)
    repository_commit = _repository_commit(REPO_ROOT)
    repository_state_files = _repository_state_files(REPO_ROOT)
    repository_state_checks = _repository_state_checks(repository_state_files)
    segment_smax = settings.target_s - settings.start_s
    tag = (
        f"{settings.nucleus}_Nrefmax{settings.nrefmax}_emax{settings.emax}"
        f"_s{path_token(settings.start_s)}to{path_token(settings.target_s)}"
        f"_tol{settings.ode_tolerance:.0e}".replace("+", "").replace("-", "m")
    )
    if settings.label:
        tag += f"_{settings.label}"
    result_dir = result_root / tag
    if result_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing MR flow: {result_dir}")
    result_dir.mkdir(parents=True)
    script_file = result_dir / f"run_{tag}.sh"
    flow_file = result_dir / f"flow_{tag}.dat"
    prefix = result_dir / tag
    output_j64 = result_dir / f"H_{tag}.jcoupled64"
    output_no2bpack = result_dir / f"H_{tag}.no2bpack"
    resource_usage = result_dir / "resource_usage.txt"
    manifest = result_dir / "metadata.json"
    nucleus_a = int(re.search(r"\d+", settings.nucleus).group())
    parameters = [
        f"2bme={interaction}",
        "fmt2=jcoupled64",
        "3bme=none",
        f"reference={settings.nucleus}",
        f"valence_space={settings.nucleus}",
        f"A={nucleus_a}",
        "hw=20",
        f"emax={settings.emax}",
        "basis=oscillator",
        f"method={settings.method}",
        "nsteps=1",
        "core_generator=white-ncsm",
        "denominator_partitioning=Epstein_Nesbet",
        "denominator_delta=0",
        "BetaCM=0",
        "nucleon_mass_correction=false",
        f"mr_reference_file={reference_file}",
        "mr_validation_tolerance=1e-10",
        f"smax={segment_smax:.17g}",
        f"ds_0={settings.ds_0:.17g}",
        f"dsmax={settings.dsmax:.17g}",
        f"ode_tolerance={settings.ode_tolerance:.17g}",
        f"eta_criterion={settings.eta_criterion:.17g}",
        f"flowfile={flow_file}",
        f"intfile={prefix}",
        f"write_H_jcoupled64={output_j64}",
        f"write_H_no2bpack={output_no2bpack}",
    ]
    command = shlex.join([str(executable_path), *parameters])
    nodelist_line = (
        f"#SBATCH --nodelist={settings.nodelist}\n" if settings.nodelist else ""
    )
    library_paths = f"{executable_path.parent}:{executable_path.parent.parent}"
    contents = f"""#!/bin/bash -l
#SBATCH --partition={settings.partition}
#SBATCH --qos={settings.qos}
#SBATCH -J mr_{settings.nucleus.lower()}_e{settings.emax}
{nodelist_line}#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task={settings.cpus}
#SBATCH -o {result_dir}/log_{tag}_%j.txt

set -eo pipefail
cd {REPO_ROOT}
source ./sourceme.sh
export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-{settings.cpus}}}
export OPENBLAS_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-{settings.cpus}}}
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:{library_paths}"

echo repository_commit={repository_commit}
{repository_state_checks}
echo cumulative_start_s={settings.start_s:.17g}
echo cumulative_target_s={settings.target_s:.17g}
echo '{interaction_sha256}  {interaction}' | sha256sum -c -
echo '{reference_sha256}  {reference_file}' | sha256sum -c -
echo '{executable_sha256}  {executable_path}' | sha256sum -c -
echo '{shared_library_sha256}  {shared_library}' | sha256sum -c -
echo '{environment_script_sha256}  {environment_script}' | sha256sum -c -
if ldd {shlex.quote(str(executable_path))} | grep 'not found'; then
  exit 1
fi
start_seconds=$SECONDS
maximum_rss_kib=0
{command} &
flow_pid=$!
while kill -0 "$flow_pid" 2>/dev/null; do
  current_rss_kib=$(awk '/^VmHWM:/ {{print $2}}' "/proc/$flow_pid/status" 2>/dev/null || true)
  if [[ "$current_rss_kib" =~ ^[0-9]+$ ]] && (( current_rss_kib > maximum_rss_kib )); then
    maximum_rss_kib=$current_rss_kib
  fi
  sleep 1
done
set +e
wait "$flow_pid"
flow_status=$?
set -e
{{
  echo wall_seconds=$((SECONDS - start_seconds))
  echo maximum_rss_kib="$maximum_rss_kib"
  echo exit_status="$flow_status"
}} | tee {shlex.quote(str(resource_usage))}
if (( flow_status != 0 )); then
  exit "$flow_status"
fi
test -s {shlex.quote(str(output_j64))}
test -s {shlex.quote(str(output_no2bpack))}
sha256sum {shlex.quote(str(output_j64))} {shlex.quote(str(output_no2bpack))}
"""
    script_file.write_text(contents, encoding="utf-8")
    script_file.chmod(0o755)
    metadata = {
        "schema": "mrimsrg_cpp_jscheme_slurm_v1",
        "repository_commit": repository_commit,
        "repository_state_files": [
            {"path": str(path), "sha256": digest}
            for path, digest in repository_state_files
        ],
        "settings": {
            **asdict(settings),
            "interaction": str(interaction),
            "reference_file": str(reference_file),
            "result_root": str(result_root),
            "executable": str(executable_path),
        },
        "cumulative_start_s": settings.start_s,
        "cumulative_target_s": settings.target_s,
        "segment_smax": segment_smax,
        "interaction_sha256": interaction_sha256,
        "reference_sha256": reference_sha256,
        "executable_sha256": executable_sha256,
        "shared_library": str(shared_library),
        "shared_library_sha256": shared_library_sha256,
        "environment_script": str(environment_script),
        "environment_script_sha256": environment_script_sha256,
        "flow_file": str(flow_file),
        "output_jcoupled64": str(output_j64),
        "output_no2bpack": str(output_no2bpack),
        "resource_usage": str(resource_usage),
        "script": str(script_file),
    }
    manifest.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    return script_file


def _validate_mr_jscheme_input_settings(settings: MRJschemeInputSettings) -> None:
    allowed_nrefmax = {"He4": (0, 2), "Be8": (0,), "C12": (0,), "O16": (0, 2)}
    if settings.nucleus not in allowed_nrefmax:
        raise ValueError("unsupported MR J-scheme input nucleus")
    if settings.nrefmax not in allowed_nrefmax[settings.nucleus]:
        raise ValueError("unsupported nucleus/Nrefmax MR input reference")
    if settings.emax < 2 or settings.emax % 2:
        raise ValueError("MR J-scheme input emax must be an even integer >= 2")
    if settings.partition not in ("c128m1024", "c128m512", "compute_C", "compute_A"):
        raise ValueError("unsupported point7 partition")
    if settings.qos != "low":
        raise ValueError("MR J-scheme input preparation requires qos=low")
    if settings.cpus < 1:
        raise ValueError("MR J-scheme input preparation cpus must be positive")
    for value, field in ((settings.nodelist, "nodelist"), (settings.label, "label")):
        if value is not None and not re.fullmatch(r"[A-Za-z0-9_-]+", value):
            raise ValueError(f"{field} contains unsupported characters")
    for path, description in (
        (settings.minipack, "minipack"),
        (settings.source_reference, "source MR reference"),
        (settings.converter, "minipack-to-J64 converter"),
        (settings.pyimsrg_dir / "pyIMSRG.so", "pyIMSRG module"),
        (REPO_ROOT / "prototype" / "mrimsrg" / "embed_jref.py", "jref embedder"),
        (REPO_ROOT / "prototype" / "mrimsrg" / "jref_format.py", "jref format module"),
        (REPO_ROOT / "prototype" / "mrimsrg" / "pyimsrg_utils.py", "pyIMSRG loader"),
        (REPO_ROOT / "sourceme.sh", "repository environment"),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"{description} is unavailable: {path}")


def generate_mr_jscheme_input_slurm(settings: MRJschemeInputSettings) -> Path:
    """Generate one hashed J64/reference preparation job for a larger MR space."""
    _validate_mr_jscheme_input_settings(settings)
    minipack = settings.minipack.resolve()
    source_reference = settings.source_reference.resolve()
    converter = settings.converter.resolve()
    pyimsrg_module = (settings.pyimsrg_dir / "pyIMSRG.so").resolve()
    pyimsrg_dir = pyimsrg_module.parent
    embedder = (REPO_ROOT / "prototype" / "mrimsrg" / "embed_jref.py").resolve()
    format_module = (REPO_ROOT / "prototype" / "mrimsrg" / "jref_format.py").resolve()
    pyimsrg_loader = (REPO_ROOT / "prototype" / "mrimsrg" / "pyimsrg_utils.py").resolve()
    environment_script = (REPO_ROOT / "sourceme.sh").resolve()
    result_root = settings.result_root.resolve()
    repository_commit = _repository_commit(REPO_ROOT)
    repository_state_files = _repository_state_files(REPO_ROOT)
    repository_state_checks = _repository_state_checks(repository_state_files)
    tag = (
        f"{settings.nucleus}_Nrefmax{settings.nrefmax}_emax{settings.emax}"
        "_input"
    )
    if settings.label:
        tag += f"_{settings.label}"
    result_dir = result_root / tag
    if result_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing MR input: {result_dir}")
    result_dir.mkdir(parents=True)
    script_file = result_dir / f"run_{tag}.sh"
    output_j64 = result_dir / f"{settings.nucleus}_emax{settings.emax}_bare.jcoupled64"
    output_reference = result_dir / f"{settings.nucleus}_emax{settings.emax}.jref"
    embed_report = result_dir / "embed_report.json"
    converter_stdout = result_dir / "converter.stdout"
    converter_time = result_dir / "converter.time"
    embed_stdout = result_dir / "embed.stdout"
    embed_time = result_dir / "embed.time"
    manifest = result_dir / "metadata.json"
    nucleus_a = int(re.search(r"\d+", settings.nucleus).group())
    hashes = {
        "minipack_sha256": _sha256(minipack),
        "source_reference_sha256": _sha256(source_reference),
        "converter_sha256": _sha256(converter),
        "pyimsrg_sha256": _sha256(pyimsrg_module),
        "embedder_sha256": _sha256(embedder),
        "format_module_sha256": _sha256(format_module),
        "pyimsrg_loader_sha256": _sha256(pyimsrg_loader),
        "environment_script_sha256": _sha256(environment_script),
    }
    nodelist_line = (
        f"#SBATCH --nodelist={settings.nodelist}\n" if settings.nodelist else ""
    )
    converter_command = shlex.join([
        str(converter), "--interaction", str(minipack),
        "--output", str(output_j64), "--A", str(nucleus_a),
    ])
    embed_command = shlex.join([
        "python3", str(embedder), "--source", str(source_reference),
        "--interaction", str(minipack), "--output", str(output_reference),
        "--pyimsrg-dir", str(pyimsrg_dir), "--json", str(embed_report),
    ])
    checks = "\n".join(
        f"echo '{digest}  {path}' | sha256sum -c -"
        for digest, path in (
            (hashes["minipack_sha256"], minipack),
            (hashes["source_reference_sha256"], source_reference),
            (hashes["converter_sha256"], converter),
            (hashes["pyimsrg_sha256"], pyimsrg_module),
            (hashes["embedder_sha256"], embedder),
            (hashes["format_module_sha256"], format_module),
            (hashes["pyimsrg_loader_sha256"], pyimsrg_loader),
            (hashes["environment_script_sha256"], environment_script),
        )
    )
    contents = f"""#!/bin/bash -l
#SBATCH --partition={settings.partition}
#SBATCH --qos={settings.qos}
#SBATCH -J mr_input_{settings.nucleus.lower()}_e{settings.emax}
{nodelist_line}#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task={settings.cpus}
#SBATCH -o {result_dir}/log_{tag}_%j.txt

set -eo pipefail
cd {REPO_ROOT}
source ./sourceme.sh
export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-{settings.cpus}}}
export OPENBLAS_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-{settings.cpus}}}
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:{pyimsrg_dir}:{converter.parent}"

echo repository_commit={repository_commit}
{repository_state_checks}
{checks}
if ldd {shlex.quote(str(converter))} | grep 'not found'; then
  exit 1
fi
run_with_rss() {{
  local command_label="$1"
  local stdout_file="$2"
  local resource_file="$3"
  shift 3
  local start_seconds=$SECONDS
  local maximum_rss_kib=0
  : > "$resource_file"
  "$@" > "$stdout_file" 2>> "$resource_file" &
  local sampled_pid=$!
  while kill -0 "$sampled_pid" 2>/dev/null; do
    local current_rss_kib
    current_rss_kib=$(awk '/^VmHWM:/ {{print $2}}' "/proc/$sampled_pid/status" 2>/dev/null || true)
    if [[ "$current_rss_kib" =~ ^[0-9]+$ ]] && (( current_rss_kib > maximum_rss_kib )); then
      maximum_rss_kib=$current_rss_kib
    fi
    sleep 1
  done
  set +e
  wait "$sampled_pid"
  local command_status=$?
  set -e
  {{
    echo command_label="$command_label"
    echo wall_seconds=$((SECONDS - start_seconds))
    echo maximum_rss_kib="$maximum_rss_kib"
    echo exit_status="$command_status"
  }} >> "$resource_file"
  return "$command_status"
}}
run_with_rss converter {shlex.quote(str(converter_stdout))} {shlex.quote(str(converter_time))} {converter_command}
run_with_rss embed_reference {shlex.quote(str(embed_stdout))} {shlex.quote(str(embed_time))} {embed_command}
test -s {shlex.quote(str(output_j64))}
test -s {shlex.quote(str(output_reference))}
test -s {shlex.quote(str(embed_report))}
sha256sum {shlex.quote(str(output_j64))} {shlex.quote(str(output_reference))} {shlex.quote(str(embed_report))}
"""
    script_file.write_text(contents, encoding="utf-8")
    script_file.chmod(0o755)
    metadata = {
        "schema": "mrimsrg_jscheme_input_slurm_v1",
        "repository_commit": repository_commit,
        "repository_state_files": [
            {"path": str(path), "sha256": digest}
            for path, digest in repository_state_files
        ],
        "settings": {
            **asdict(settings),
            "minipack": str(minipack),
            "source_reference": str(source_reference),
            "result_root": str(result_root),
            "converter": str(converter),
            "pyimsrg_dir": str(pyimsrg_dir),
        },
        **hashes,
        "pyimsrg_module": str(pyimsrg_module),
        "embedder": str(embedder),
        "format_module": str(format_module),
        "pyimsrg_loader": str(pyimsrg_loader),
        "environment_script": str(environment_script),
        "output_jcoupled64": str(output_j64),
        "output_reference": str(output_reference),
        "embed_report": str(embed_report),
        "converter_stdout": str(converter_stdout),
        "converter_time": str(converter_time),
        "embed_stdout": str(embed_stdout),
        "embed_time": str(embed_time),
        "script": str(script_file),
    }
    manifest.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    return script_file


def generate_mr_ncsm_readback_slurm(settings: MRNCSMReadbackSettings) -> Path:
    """Generate one hashed downstream NCSM/no2bpack readback job."""
    if min(settings.proton_number, settings.neutron_number, settings.nmax) < 0:
        raise ValueError("MR NCSM particle numbers and Nmax must be nonnegative")
    if settings.proton_number + settings.neutron_number < 2:
        raise ValueError("MR NCSM readback requires at least two particles")
    if settings.nmax % 2:
        raise ValueError("MR NCSM readback currently requires even Nmax")
    if min(settings.states, settings.max_iter, settings.cpus) < 1:
        raise ValueError("MR NCSM states, iterations, and cpus must be positive")
    if settings.partition not in ("c128m1024", "c128m512", "compute_C", "compute_A"):
        raise ValueError("unsupported point7 partition")
    if settings.qos != "low":
        raise ValueError("MR NCSM production requires qos=low")
    for value, field in ((settings.nodelist, "nodelist"), (settings.label, "label")):
        if value is not None and not re.fullmatch(r"[A-Za-z0-9_-]+", value):
            raise ValueError(f"{field} contains unsupported characters")
    if not settings.no2bpack.is_file():
        raise FileNotFoundError(f"no2bpack is unavailable: {settings.no2bpack}")
    if not settings.executable.is_file():
        raise FileNotFoundError(f"NCSM validator is unavailable: {settings.executable}")
    environment_script = (REPO_ROOT / "sourceme.sh").resolve()
    if not environment_script.is_file():
        raise FileNotFoundError("repository sourceme.sh is unavailable")

    no2bpack = settings.no2bpack.resolve()
    executable = settings.executable.resolve()
    result_root = settings.result_root.resolve()
    repository_commit = _repository_commit(REPO_ROOT)
    repository_state_files = _repository_state_files(REPO_ROOT)
    repository_state_checks = _repository_state_checks(repository_state_files)
    no2bpack_sha256 = _sha256(no2bpack)
    executable_sha256 = _sha256(executable)
    environment_script_sha256 = _sha256(environment_script)
    tag = (
        f"Z{settings.proton_number}_N{settings.neutron_number}"
        f"_Nmax{settings.nmax}_states{settings.states}"
    )
    if settings.label:
        tag += f"_{settings.label}"
    result_dir = result_root / tag
    if result_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing MR NCSM readback: {result_dir}")
    result_dir.mkdir(parents=True)
    script_file = result_dir / f"run_{tag}.sh"
    log_file = result_dir / f"log_{tag}_%j.txt"
    manifest = result_dir / "metadata.json"
    command = shlex.join([
        str(executable),
        "--no2bpack", str(no2bpack),
        "--Z", str(settings.proton_number),
        "--N", str(settings.neutron_number),
        "--nmax", str(settings.nmax),
        "--states", str(settings.states),
        "--max-iter", str(settings.max_iter),
    ])
    if settings.sample_rss:
        execution_block = f"""{command} &
ncsm_pid=$!
maximum_rss_kib=0
while kill -0 \"$ncsm_pid\" 2>/dev/null; do
  current_rss_kib=$(awk '/^VmHWM:/ {{print $2}}' \"/proc/$ncsm_pid/status\" 2>/dev/null || true)
  if [[ \"$current_rss_kib\" =~ ^[0-9]+$ ]] && (( current_rss_kib > maximum_rss_kib )); then
    maximum_rss_kib=$current_rss_kib
  fi
  sleep 1
done
set +e
wait \"$ncsm_pid\"
ncsm_status=$?
set -e
echo wall_seconds=$((SECONDS - start_seconds))
echo maximum_rss_kib=$maximum_rss_kib
exit \"$ncsm_status\""""
    else:
        execution_block = f"""set +e
{command}
ncsm_status=$?
set -e
echo wall_seconds=$((SECONDS - start_seconds))
echo maximum_rss_kib=not_sampled
exit \"$ncsm_status\""""
    nodelist_line = (
        f"#SBATCH --nodelist={settings.nodelist}\n" if settings.nodelist else ""
    )
    contents = f"""#!/bin/bash -l
#SBATCH --partition={settings.partition}
#SBATCH --qos={settings.qos}
#SBATCH -J mr_ncsm_A{settings.proton_number + settings.neutron_number}
{nodelist_line}#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task={settings.cpus}
#SBATCH -o {log_file}

set -eo pipefail
cd {REPO_ROOT}
source ./sourceme.sh
export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-{settings.cpus}}}
export OPENBLAS_NUM_THREADS=1

echo repository_commit={repository_commit}
{repository_state_checks}
echo '{no2bpack_sha256}  {no2bpack}' | sha256sum -c -
echo '{executable_sha256}  {executable}' | sha256sum -c -
echo '{environment_script_sha256}  {environment_script}' | sha256sum -c -
if ldd {shlex.quote(str(executable))} | grep 'not found'; then
  exit 1
fi
start_seconds=$SECONDS
{execution_block}
"""
    script_file.write_text(contents, encoding="utf-8")
    script_file.chmod(0o755)
    metadata = {
        "schema": "mrimsrg_ncsm_readback_slurm_v1",
        "repository_commit": repository_commit,
        "repository_state_files": [
            {"path": str(path), "sha256": digest}
            for path, digest in repository_state_files
        ],
        "settings": {
            **asdict(settings),
            "no2bpack": str(no2bpack),
            "result_root": str(result_root),
            "executable": str(executable),
        },
        "no2bpack_sha256": no2bpack_sha256,
        "executable_sha256": executable_sha256,
        "environment_script": str(environment_script),
        "environment_script_sha256": environment_script_sha256,
        "script": str(script_file),
        "log_pattern": str(log_file),
    }
    manifest.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    return script_file


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
    parser.add_argument(
        "--mr-jscheme", action="store_true",
        help="generate one production C++ J-scheme MR direct-flow job",
    )
    parser.add_argument(
        "--mr-prepare-jscheme", action="store_true",
        help="prepare one larger-space J64 interaction and embedded MR reference",
    )
    parser.add_argument(
        "--mr-ncsm-readback", action="store_true",
        help="generate one downstream NCSM/no2bpack readback job",
    )
    parser.add_argument("--mr-nucleus", choices=("He4", "Be8", "C12", "O16"))
    parser.add_argument("--mr-nrefmax", type=int, choices=(0, 2))
    parser.add_argument("--mr-emax", type=int)
    parser.add_argument("--mr-interaction", type=Path)
    parser.add_argument("--mr-reference-file", type=Path)
    parser.add_argument("--mr-minipack", type=Path)
    parser.add_argument("--mr-source-reference", type=Path)
    parser.add_argument("--mr-converter", type=Path)
    parser.add_argument("--mr-pyimsrg-dir", type=Path)
    parser.add_argument("--mr-start-s", type=float, default=0.0)
    parser.add_argument("--mr-target-s", type=float)
    parser.add_argument(
        "--mr-method", choices=("flow", "flow_adaptive", "flow_RK4"),
        default="flow",
    )
    parser.add_argument("--mr-ds-0", type=float, default=1e-2)
    parser.add_argument("--mr-dsmax", type=float, default=0.5)
    parser.add_argument("--mr-ode-tolerance", type=float, default=1e-8)
    parser.add_argument("--mr-eta-criterion", type=float, default=1e-6)
    parser.add_argument(
        "--mr-partition",
        choices=("c128m1024", "c128m512", "compute_C", "compute_A"),
        default="c128m512",
    )
    parser.add_argument("--mr-nodelist")
    parser.add_argument("--mr-cpus", type=int, default=64)
    parser.add_argument("--mr-label")
    parser.add_argument("--mr-result-root", type=Path)
    parser.add_argument("--mr-executable", type=Path)
    parser.add_argument("--mr-ncsm-no2bpack", type=Path)
    parser.add_argument("--mr-z", type=int)
    parser.add_argument("--mr-n", type=int)
    parser.add_argument("--mr-nmax", type=int)
    parser.add_argument("--mr-states", type=int, default=3)
    parser.add_argument("--mr-max-iter", type=int, default=300)
    parser.add_argument(
        "--mr-ncsm-no-rss-sampling", action="store_true",
        help="run the NCSM validator in the foreground without /proc RSS sampling",
    )
    parser.add_argument("--mr-ncsm-executable", type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.generate_only and args.submit:
        raise SystemExit("--generate-only and --submit are mutually exclusive")
    mr_modes = sum((args.mr_jscheme, args.mr_prepare_jscheme,
                    args.mr_ncsm_readback))
    if mr_modes > 1:
        raise SystemExit("MR J-scheme modes are mutually exclusive")
    check_and_make_dir(REPO_ROOT / "result")

    if args.smoke_test:
        run_smoke_test()
        exit(0)

    if args.mr_jscheme:
        required = {
            "--mr-nucleus": args.mr_nucleus,
            "--mr-nrefmax": args.mr_nrefmax,
            "--mr-emax": args.mr_emax,
            "--mr-interaction": args.mr_interaction,
            "--mr-reference-file": args.mr_reference_file,
            "--mr-target-s": args.mr_target_s,
        }
        missing = [name for name, value in required.items() if value is None]
        if missing:
            raise SystemExit("--mr-jscheme requires " + ", ".join(missing))
        default_executable = REPO_ROOT / "build" / "src" / "imsrg++"
        if not default_executable.is_file():
            default_executable = REPO_ROOT / "build" / "imsrg++"
        script_file = generate_mr_jscheme_slurm(
            MRJschemeSettings(
                nucleus=args.mr_nucleus,
                nrefmax=args.mr_nrefmax,
                emax=args.mr_emax,
                interaction=args.mr_interaction,
                reference_file=args.mr_reference_file,
                start_s=args.mr_start_s,
                target_s=args.mr_target_s,
                method=args.mr_method,
                ds_0=args.mr_ds_0,
                dsmax=args.mr_dsmax,
                ode_tolerance=args.mr_ode_tolerance,
                eta_criterion=args.mr_eta_criterion,
                partition=args.mr_partition,
                cpus=args.mr_cpus,
                nodelist=args.mr_nodelist,
                label=args.mr_label,
                result_root=(
                    args.mr_result_root
                    if args.mr_result_root is not None
                    else REPO_ROOT / "result" / "mr-jscheme-flow"
                ),
                executable=(
                    args.mr_executable
                    if args.mr_executable is not None
                    else default_executable
                ),
            )
        )
        print(script_file)
        if args.generate_only:
            exit(0)
        if args.submit:
            result = subprocess.run(
                ["sbatch", str(script_file)], capture_output=True, text=True, check=True
            )
            if result.stdout.strip():
                print(result.stdout.strip())
            exit(0)
        run_mode = input("submit job: y/n\n")
        if run_mode.lower()[0] == "y":
            subprocess.run(["sbatch", str(script_file)], check=True)
        else:
            subprocess.run(["bash", str(script_file)], check=True)
        exit(0)

    if args.mr_prepare_jscheme:
        required = {
            "--mr-nucleus": args.mr_nucleus,
            "--mr-nrefmax": args.mr_nrefmax,
            "--mr-emax": args.mr_emax,
            "--mr-minipack": args.mr_minipack,
            "--mr-source-reference": args.mr_source_reference,
        }
        missing = [name for name, value in required.items() if value is None]
        if missing:
            raise SystemExit("--mr-prepare-jscheme requires " + ", ".join(missing))
        script_file = generate_mr_jscheme_input_slurm(
            MRJschemeInputSettings(
                nucleus=args.mr_nucleus,
                nrefmax=args.mr_nrefmax,
                emax=args.mr_emax,
                minipack=args.mr_minipack,
                source_reference=args.mr_source_reference,
                partition=args.mr_partition,
                cpus=args.mr_cpus,
                nodelist=args.mr_nodelist,
                label=args.mr_label,
                result_root=(
                    args.mr_result_root
                    if args.mr_result_root is not None
                    else REPO_ROOT / "result" / "mr-jscheme-inputs"
                ),
                converter=(
                    args.mr_converter
                    if args.mr_converter is not None
                    else REPO_ROOT / "prototype" / "mrimsrg" / "build" /
                    "mrimsrg_minipack_to_j64"
                ),
                pyimsrg_dir=(
                    args.mr_pyimsrg_dir
                    if args.mr_pyimsrg_dir is not None
                    else REPO_ROOT / "build" / "src"
                ),
            )
        )
        print(script_file)
        if args.generate_only:
            exit(0)
        if args.submit:
            result = subprocess.run(
                ["sbatch", str(script_file)], capture_output=True, text=True, check=True
            )
            if result.stdout.strip():
                print(result.stdout.strip())
            exit(0)
        run_mode = input("submit job: y/n\n")
        if run_mode.lower()[0] == "y":
            subprocess.run(["sbatch", str(script_file)], check=True)
        else:
            subprocess.run(["bash", str(script_file)], check=True)
        exit(0)

    if args.mr_ncsm_readback:
        required = {
            "--mr-ncsm-no2bpack": args.mr_ncsm_no2bpack,
            "--mr-z": args.mr_z,
            "--mr-n": args.mr_n,
            "--mr-nmax": args.mr_nmax,
        }
        missing = [name for name, value in required.items() if value is None]
        if missing:
            raise SystemExit("--mr-ncsm-readback requires " + ", ".join(missing))
        script_file = generate_mr_ncsm_readback_slurm(
            MRNCSMReadbackSettings(
                no2bpack=args.mr_ncsm_no2bpack,
                proton_number=args.mr_z,
                neutron_number=args.mr_n,
                nmax=args.mr_nmax,
                states=args.mr_states,
                max_iter=args.mr_max_iter,
                partition=args.mr_partition,
                cpus=args.mr_cpus,
                nodelist=args.mr_nodelist,
                label=args.mr_label,
                sample_rss=not args.mr_ncsm_no_rss_sampling,
                result_root=(
                    args.mr_result_root
                    if args.mr_result_root is not None
                    else REPO_ROOT / "result" / "mr-ncsm-readback"
                ),
                executable=(
                    args.mr_ncsm_executable
                    if args.mr_ncsm_executable is not None
                    else REPO_ROOT / "prototype" / "mrimsrg" / "build-emax6" / "mrimsrg_validate"
                ),
            )
        )
        print(script_file)
        if args.generate_only:
            exit(0)
        if args.submit:
            result = subprocess.run(
                ["sbatch", str(script_file)], capture_output=True, text=True, check=True
            )
            if result.stdout.strip():
                print(result.stdout.strip())
            exit(0)
        run_mode = input("submit job: y/n\n")
        if run_mode.lower()[0] == "y":
            subprocess.run(["sbatch", str(script_file)], check=True)
        else:
            subprocess.run(["bash", str(script_file)], check=True)
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
