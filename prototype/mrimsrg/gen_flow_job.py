"""Generate one inspected point7 Slurm job for the rapid MR-IMSRG flow.

This is intentionally a single-parameter-set generator.  Invoke it once for
each nucleus or tolerance; it never creates or submits a scan.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import re
import subprocess


_REFERENCE_NAMES = {
    ("He4", 0): "He4_Nrefmax0_final",
    ("Be8", 0): "Be8_Nrefmax0_final",
    ("C12", 0): "C12_Nrefmax0_final",
    ("O16", 0): "O16_Nrefmax0_final",
    ("He4", 2): "He4_Nrefmax2",
    ("O16", 2): "O16_Nrefmax2",
}
_PARTITIONS = ("c128m1024", "c128m512", "compute_C", "compute_A")
_POINT7_INTERACTION = Path(
    "/tns/mengziyan/mr-imsrg-inputs/"
    "TwBME_N2LO_opt_hw20_emax2_e2max4.minipack"
)
_POINT7_PYTHON = Path("/opt/library/miniconda-3.12.9/bin/python3")
_POINT7_SYSTEM_PYTHON = Path("/usr/bin/python3")


@dataclass(frozen=True)
class JobSettings:
    nucleus: str
    nrefmax: int = 0
    interaction: Path = _POINT7_INTERACTION
    label: str | None = None
    partition: str = "c128m512"
    nodelist: str | None = None
    smax: float = 25000.0
    checkpoint_s: float = 40.0
    rtol: float = 1e-6
    atol: float = 1e-8
    max_step: float = 10.0
    residual_ratio: float = 1e-6
    resume_from: Path | None = None


@dataclass(frozen=True)
class SRCheckSettings:
    nucleus: str
    production_flow: Path
    jcoupled64: Path
    pyimsrg_dir: Path
    nrefmax: int = 0
    interaction: Path = _POINT7_INTERACTION
    label: str | None = None
    partition: str = "c128m512"
    nodelist: str | None = None
    ode_tolerance: float = 1e-8
    initial_step: float = 1e-2
    flow_tolerance: float = 1e-5


@dataclass(frozen=True)
class MRCheckSettings:
    nucleus: str
    nrefmax: int
    production_flow: Path
    pyimsrg_dir: Path
    interaction: Path = _POINT7_INTERACTION
    label: str | None = None
    partition: str = "c128m512"
    nodelist: str | None = None
    ode_tolerance: float = 1e-9
    initial_step: float = 1e-2
    flow_tolerance: float = 1e-5
    python_executable: Path = _POINT7_SYSTEM_PYTHON


def _number_tag(value: float) -> str:
    return format(value, ".0e").replace("+", "").replace("-", "m")


def _repository_commit(repo_root: Path) -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        return "unknown"
    return completed.stdout.strip()


def _validate(settings: JobSettings) -> None:
    if (settings.nucleus, settings.nrefmax) not in _REFERENCE_NAMES:
        raise ValueError(
            f"unsupported reference: {settings.nucleus}, Nrefmax={settings.nrefmax}"
        )
    if settings.partition not in _PARTITIONS:
        raise ValueError(f"unsupported point7 partition: {settings.partition}")
    if settings.nodelist is not None and not re.fullmatch(r"[A-Za-z0-9_-]+", settings.nodelist):
        raise ValueError("nodelist contains unsupported characters")
    if settings.label is not None and not re.fullmatch(r"[A-Za-z0-9_-]+", settings.label):
        raise ValueError("label contains unsupported characters")
    start_s = 0.0
    if settings.resume_from is not None:
        metadata_path = settings.resume_from / "metadata.json"
        if not metadata_path.is_file():
            raise FileNotFoundError(
                f"resume flow is unavailable: {settings.resume_from}"
            )
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        trajectory = metadata.get("trajectory")
        if not isinstance(trajectory, list) or not trajectory:
            raise ValueError("resume flow has no recorded trajectory")
        start_s = float(trajectory[-1]["s"])
    if not start_s < settings.checkpoint_s < settings.smax:
        raise ValueError("checkpoint_s must lie strictly between start_s and smax")
    if settings.rtol <= 0.0 or settings.atol <= 0.0 or settings.max_step <= 0.0:
        raise ValueError("ODE tolerances and max_step must be positive")
    if not 0.0 < settings.residual_ratio < 1.0:
        raise ValueError("residual_ratio must lie between zero and one")
    if not settings.interaction.is_file():
        raise FileNotFoundError(f"fixed interaction is unavailable: {settings.interaction}")


def generate_job(
    repo_root: Path, result_root: Path, settings: JobSettings
) -> Path:
    """Generate exactly one job script and return its absolute path."""
    _validate(settings)
    repo_root = repo_root.resolve()
    result_root = result_root.resolve()
    repository_commit = _repository_commit(repo_root)
    reference = (
        repo_root
        / "prototype"
        / "mrimsrg"
        / "data"
        / _REFERENCE_NAMES[(settings.nucleus, settings.nrefmax)]
    )
    if not (reference / "metadata.json").is_file():
        raise FileNotFoundError(f"reference bundle is unavailable: {reference}")

    tag = (
        f"{settings.nucleus}_Nrefmax{settings.nrefmax}_rtol{_number_tag(settings.rtol)}"
        f"_atol{_number_tag(settings.atol)}"
    )
    if settings.label:
        tag += f"_{settings.label}"
    case_root = result_root / tag
    output = case_root / "flow"
    if output.exists():
        raise FileExistsError(f"refusing to overwrite an existing flow: {output}")
    case_root.mkdir(parents=True, exist_ok=True)
    script = case_root / f"job_{tag}.sh"

    nodelist_line = (
        f"#SBATCH --nodelist={settings.nodelist}\n" if settings.nodelist else ""
    )
    command = (
        f"PYTHONPATH={repo_root / 'prototype' / 'mrimsrg'} {_POINT7_PYTHON} -u "
        f"{repo_root / 'prototype' / 'mrimsrg' / 'run_mrimsrg.py'} "
        f"{reference} {output} "
        f"--interaction {settings.interaction.resolve()} "
        f"--smax {settings.smax:.17g} "
        f"--checkpoint-s {settings.checkpoint_s:.17g} "
        f"--rtol {settings.rtol:.17g} --atol {settings.atol:.17g} "
        f"--max-step {settings.max_step:.17g} "
        f"--residual-ratio {settings.residual_ratio:.17g}"
    )
    if settings.resume_from is not None:
        command += f" --resume-from {settings.resume_from.resolve()}"
    contents = f"""#!/bin/bash
#SBATCH --job-name=mr_{settings.nucleus.lower()}
#SBATCH --partition={settings.partition}
#SBATCH --qos=low
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
{nodelist_line}#SBATCH -o {case_root}/log_{tag}_%j.txt

set -eo pipefail
cd {repo_root}
source /opt/modules/init/bash
source ./sourceme.sh
export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-64}}
export OPENBLAS_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-64}}
export MKL_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-64}}

{_POINT7_PYTHON} -c 'import sys, numpy, scipy; print("python", sys.version.split()[0], "numpy", numpy.__version__, "scipy", scipy.__version__, flush=True)'
echo repository_commit={repository_commit}
{command}
"""
    script.write_text(contents, encoding="utf-8")
    script.chmod(0o755)
    return script


def generate_sr_check_job(
    repo_root: Path, result_root: Path, settings: SRCheckSettings
) -> Path:
    """Generate one full-flow SR degeneration comparison job."""
    common = JobSettings(
        nucleus=settings.nucleus,
        nrefmax=settings.nrefmax,
        interaction=settings.interaction,
        label=settings.label,
        partition=settings.partition,
        nodelist=settings.nodelist,
    )
    _validate(common)
    if settings.nrefmax != 0 or settings.nucleus not in ("He4", "O16"):
        raise ValueError("strict full-flow SR checks support He4/O16 Nrefmax=0")
    if not (settings.production_flow / "metadata.json").is_file():
        raise FileNotFoundError(
            f"production flow is unavailable: {settings.production_flow}"
        )
    if not settings.jcoupled64.is_file():
        raise FileNotFoundError(f"jcoupled64 is unavailable: {settings.jcoupled64}")
    if not (settings.pyimsrg_dir / "pyIMSRG.so").is_file():
        raise FileNotFoundError(f"pyIMSRG module is unavailable: {settings.pyimsrg_dir}")
    if min(settings.ode_tolerance, settings.initial_step, settings.flow_tolerance) <= 0:
        raise ValueError("SR check tolerances and initial step must be positive")

    repo_root = repo_root.resolve()
    result_root = result_root.resolve()
    repository_commit = _repository_commit(repo_root)
    reference = (
        repo_root
        / "prototype"
        / "mrimsrg"
        / "data"
        / _REFERENCE_NAMES[(settings.nucleus, settings.nrefmax)]
    )
    if not (reference / "metadata.json").is_file():
        raise FileNotFoundError(f"reference bundle is unavailable: {reference}")
    tag = (
        f"{settings.nucleus}_Nrefmax0_srcheck_"
        f"tol{_number_tag(settings.ode_tolerance)}"
    )
    if settings.label:
        tag += f"_{settings.label}"
    case_root = result_root / tag
    report = case_root / "report.json"
    if case_root.exists():
        raise FileExistsError(f"refusing to overwrite an existing check: {case_root}")
    case_root.mkdir(parents=True)
    script = case_root / f"job_{tag}.sh"
    nodelist_line = (
        f"#SBATCH --nodelist={settings.nodelist}\n" if settings.nodelist else ""
    )
    command = (
        f"PYTHONPATH={repo_root / 'prototype' / 'mrimsrg'} {_POINT7_PYTHON} -u "
        f"{repo_root / 'prototype' / 'mrimsrg' / 'sr_imsrgpp_check.py'} "
        f"--nucleus {settings.nucleus} --reference {reference} "
        f"--jcoupled64 {settings.jcoupled64.resolve()} "
        f"--interaction {settings.interaction.resolve()} "
        f"--pyimsrg-dir {settings.pyimsrg_dir.resolve()} "
        f"--imsrg-source {repo_root} "
        f"--production-flow {settings.production_flow.resolve()} "
        f"--full-flow-ode-tolerance {settings.ode_tolerance:.17g} "
        f"--full-flow-initial-step {settings.initial_step:.17g} "
        f"--full-flow-tolerance {settings.flow_tolerance:.17g} "
        f"--full-flow-imsrgpp-jcoupled64 {case_root / 'imsrgpp_vacuum.jcoupled64'} "
        f"--json {report}"
    )
    contents = f"""#!/bin/bash
#SBATCH --job-name=sr_{settings.nucleus.lower()}
#SBATCH --partition={settings.partition}
#SBATCH --qos=low
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
{nodelist_line}#SBATCH -o {case_root}/log_{tag}_%j.txt

set -eo pipefail
cd {repo_root}
source /opt/modules/init/bash
source ./sourceme.sh
export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-64}}
export OPENBLAS_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-64}}

echo repository_commit={repository_commit}
c++ --version | head -1
ldd {settings.pyimsrg_dir.resolve() / 'pyIMSRG.so'} | grep -E 'openblas|blas|lapack|gsl|boost|not found' || true
{_POINT7_PYTHON} -c 'import sys, numpy, scipy, sympy; print("python", sys.version.split()[0], "numpy", numpy.__version__, "scipy", scipy.__version__, "sympy", sympy.__version__, flush=True)'
{command}
"""
    script.write_text(contents, encoding="utf-8")
    script.chmod(0o755)
    return script


def generate_mr_check_job(
    repo_root: Path, result_root: Path, settings: MRCheckSettings
) -> Path:
    """Generate one correlated-reference C++/Python fixed-s comparison job."""
    common = JobSettings(
        nucleus=settings.nucleus,
        nrefmax=settings.nrefmax,
        interaction=settings.interaction,
        label=settings.label,
        partition=settings.partition,
        nodelist=settings.nodelist,
    )
    _validate(common)
    if not (settings.production_flow / "metadata.json").is_file():
        raise FileNotFoundError(
            f"production flow is unavailable: {settings.production_flow}"
        )
    if not (settings.pyimsrg_dir / "pyIMSRG.so").is_file():
        raise FileNotFoundError(f"pyIMSRG module is unavailable: {settings.pyimsrg_dir}")
    if min(settings.ode_tolerance, settings.initial_step, settings.flow_tolerance) <= 0:
        raise ValueError("MR check tolerances and initial step must be positive")
    if not settings.python_executable.is_file():
        raise FileNotFoundError(
            f"MR check Python is unavailable: {settings.python_executable}"
        )

    repo_root = repo_root.resolve()
    result_root = result_root.resolve()
    reference = (
        repo_root / "prototype" / "mrimsrg" / "data"
        / _REFERENCE_NAMES[(settings.nucleus, settings.nrefmax)]
    )
    if not (reference / "metadata.json").is_file():
        raise FileNotFoundError(f"reference bundle is unavailable: {reference}")
    tag = (
        f"{settings.nucleus}_Nrefmax{settings.nrefmax}_mrcheck_"
        f"tol{_number_tag(settings.ode_tolerance)}"
    )
    if settings.label:
        tag += f"_{settings.label}"
    case_root = result_root / tag
    if case_root.exists():
        raise FileExistsError(f"refusing to overwrite an existing check: {case_root}")
    case_root.mkdir(parents=True)
    script = case_root / f"job_{tag}.sh"
    nodelist_line = (
        f"#SBATCH --nodelist={settings.nodelist}\n" if settings.nodelist else ""
    )
    report = case_root / "report.json"
    output_j64 = case_root / "cpp_vacuum.jcoupled64"
    command = (
        f"PYTHONPATH={repo_root / 'prototype' / 'mrimsrg'} "
        f"{settings.python_executable.resolve()} -u "
        f"{repo_root / 'prototype' / 'mrimsrg' / 'mr_imsrgpp_flow_check.py'} "
        f"--nucleus {settings.nucleus} --reference {reference} "
        f"--interaction {settings.interaction.resolve()} "
        f"--production-flow {settings.production_flow.resolve()} "
        f"--pyimsrg-dir {settings.pyimsrg_dir.resolve()} "
        f"--ode-tolerance {settings.ode_tolerance:.17g} "
        f"--initial-step {settings.initial_step:.17g} "
        f"--flow-tolerance {settings.flow_tolerance:.17g} "
        f"--output-jcoupled64 {output_j64} --json {report}"
    )
    contents = f"""#!/bin/bash
#SBATCH --job-name=mrcheck_{settings.nucleus.lower()}
#SBATCH --partition={settings.partition}
#SBATCH --qos=low
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
{nodelist_line}#SBATCH -o {case_root}/log_{tag}_%j.txt

set -eo pipefail
cd {repo_root}
source /opt/modules/init/bash
source ./sourceme.sh
export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-64}}
export OPENBLAS_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-64}}

echo repository_commit={_repository_commit(repo_root)}
ldd {settings.pyimsrg_dir.resolve() / 'pyIMSRG.so'} | grep -E 'openblas|blas|lapack|gsl|boost|not found' || true
{settings.python_executable.resolve()} -c 'import sys, numpy, scipy; print("python", sys.version.split()[0], "numpy", numpy.__version__, "scipy", scipy.__version__, flush=True)'
{command}
"""
    script.write_text(contents, encoding="utf-8")
    script.chmod(0o755)
    return script


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--nucleus",
        required=True,
        choices=tuple(sorted({nucleus for nucleus, _ in _REFERENCE_NAMES})),
    )
    parser.add_argument("--nrefmax", type=int, choices=(0, 2), default=0)
    parser.add_argument("--interaction", type=Path, default=_POINT7_INTERACTION)
    parser.add_argument(
        "--label",
        help="optional single-run output label, for example smax5",
    )
    parser.add_argument("--partition", choices=_PARTITIONS, default="c128m512")
    parser.add_argument("--nodelist")
    parser.add_argument("--smax", type=float, default=25000.0)
    parser.add_argument("--checkpoint-s", type=float, default=40.0)
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--atol", type=float, default=1e-8)
    parser.add_argument("--max-step", type=float, default=10.0)
    parser.add_argument("--residual-ratio", type=float, default=1e-6)
    parser.add_argument("--resume-from", type=Path)
    parser.add_argument(
        "--sr-check-flow",
        type=Path,
        help="generate a strict C++/production comparison for this flow instead",
    )
    parser.add_argument(
        "--mr-check-flow",
        type=Path,
        help="generate a correlated-reference C++/Python fixed-s comparison",
    )
    parser.add_argument("--jcoupled64", type=Path)
    parser.add_argument("--pyimsrg-dir", type=Path)
    parser.add_argument(
        "--mr-check-python",
        type=Path,
        default=_POINT7_SYSTEM_PYTHON,
        help="Python executable ABI-matched to --pyimsrg-dir",
    )
    parser.add_argument("--sr-check-ode-tolerance", type=float)
    parser.add_argument("--sr-check-initial-step", type=float, default=1e-2)
    parser.add_argument("--sr-check-flow-tolerance", type=float, default=1e-5)
    parser.add_argument("--submit", action="store_true")
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--result-root", type=Path)
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    result_root = (
        args.result_root.resolve()
        if args.result_root is not None
        else repo_root / "result" / "mrimsrg-flow"
    )
    settings = JobSettings(
        nucleus=args.nucleus,
        nrefmax=args.nrefmax,
        interaction=args.interaction,
        label=args.label,
        partition=args.partition,
        nodelist=args.nodelist,
        smax=args.smax,
        checkpoint_s=args.checkpoint_s,
        rtol=args.rtol,
        atol=args.atol,
        max_step=args.max_step,
        residual_ratio=args.residual_ratio,
        resume_from=args.resume_from,
    )
    if args.sr_check_flow is not None and args.mr_check_flow is not None:
        parser.error("--sr-check-flow and --mr-check-flow are mutually exclusive")
    if args.sr_check_flow is None and args.mr_check_flow is None:
        script = generate_job(repo_root, result_root, settings)
    elif args.sr_check_flow is not None:
        if args.jcoupled64 is None or args.pyimsrg_dir is None:
            parser.error("--sr-check-flow requires --jcoupled64 and --pyimsrg-dir")
        ode_tolerance = (
            args.sr_check_ode_tolerance
            if args.sr_check_ode_tolerance is not None
            else args.rtol
        )
        script = generate_sr_check_job(
            repo_root,
            result_root,
            SRCheckSettings(
                nucleus=args.nucleus,
                nrefmax=args.nrefmax,
                interaction=args.interaction,
                label=args.label,
                partition=args.partition,
                nodelist=args.nodelist,
                production_flow=args.sr_check_flow,
                jcoupled64=args.jcoupled64,
                pyimsrg_dir=args.pyimsrg_dir,
                ode_tolerance=ode_tolerance,
                initial_step=args.sr_check_initial_step,
                flow_tolerance=args.sr_check_flow_tolerance,
            ),
        )
    else:
        if args.pyimsrg_dir is None:
            parser.error("--mr-check-flow requires --pyimsrg-dir")
        ode_tolerance = (
            args.sr_check_ode_tolerance
            if args.sr_check_ode_tolerance is not None
            else args.rtol
        )
        script = generate_mr_check_job(
            repo_root,
            result_root,
            MRCheckSettings(
                nucleus=args.nucleus,
                nrefmax=args.nrefmax,
                interaction=args.interaction,
                label=args.label,
                partition=args.partition,
                nodelist=args.nodelist,
                production_flow=args.mr_check_flow,
                pyimsrg_dir=args.pyimsrg_dir,
                ode_tolerance=ode_tolerance,
                initial_step=args.sr_check_initial_step,
                flow_tolerance=args.sr_check_flow_tolerance,
                python_executable=args.mr_check_python,
            ),
        )
    print(f"generated {script}")
    print(script.read_text(encoding="utf-8"), end="")
    if args.submit:
        completed = subprocess.run(
            ["sbatch", str(script)], check=True, text=True, capture_output=True
        )
        print(completed.stdout, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
