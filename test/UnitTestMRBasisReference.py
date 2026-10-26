#!/usr/bin/env python3

"""Reference computed in the imsrg++ HF basis (IM-NCSM with an HF basis).

Chain: imsrg++ single-reference basis=HF exports the HF-basis Hamiltonian
(no2bpack, 3N in NO2B with respect to HF) and the HO -> HF matrix C.mat;
BIGSTICK (no2 reader) solves the Nrefmax=0 reference in that basis;
bigstick_jref.py composes C_total = C_HF C_nat; the MR driver reads the
HO-basis NN + 3N with that reference.  For closed-shell O16 (emax2) the
reference is the HF determinant, so

1. HF energy = BIGSTICK reference energy = MR initial zero body (with the
   MR-NO2B of the 3N, which equals the HF NO2B here);
2. the MR White-NCSM flow equals imsrg++ single-reference basis=HF with the
   White generator (s = 1);
3. the export in the MR basis (mr_export_basis=natural) gives back E(s) as the
   BIGSTICK Nmax=0 energy (BIGSTICK stores float32 matrix elements: 1e-4).

The 3N file (N3LO lnl, emax6, hw16) only supplies realistic matrix elements;
the energies are not physical (hw differs from the NN part).
Needs BIGSTICK (env BIGSTICK or ~/BigstickPublick/bin/bigstick-openmp.x).
"""

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import pyIMSRG

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))
sys.path.insert(0, str(REPOSITORY / "test"))

from prototype.mrimsrg.bigstick_jref import convert  # noqa: E402
from prototype.mrimsrg.reference_io import load_reference  # noqa: E402
from UnitTestMRDriver import vacuum_operator  # noqa: E402

THREE_BODY = (REPOSITORY / "prototype/mrimsrg/data/N3LO_lnl_3N/"
              "ThBME_NO2B_srg2.0_N3LO_EM500_lnl_N3max40_hw16_emax6_e2max12_e3max6.me3j.gz")
BIGSTICK = Path(os.environ.get("BIGSTICK", Path.home() / "BigstickPublick/bin/bigstick-openmp.x"))


def run(command, log):
    environment = dict(os.environ, OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="1")
    with open(log, "w") as stream:
        status = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, env=environment).returncode
    if status != 0:
        print(Path(log).read_text()[-3000:])
        raise RuntimeError(f"command failed: {command[0]}")
    return Path(log).read_text()


def bigstick(directory, deck, log):
    with open(directory / log, "w") as stream:
        subprocess.run([str(BIGSTICK)], input=deck, text=True, stdout=stream, stderr=subprocess.STDOUT,
                       cwd=directory, check=True)


def first_energy(res):
    lines = res.read_text().splitlines()
    index = next(i for i, line in enumerate(lines) if line.split()[:3] == ["State", "E", "Ex"])
    return float(lines[index + 1].split()[1])


def last_flow_energy(path):
    rows = [line.split() for line in Path(path).read_text().splitlines() if line.strip()]
    return float(rows[-1][2])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("imsrg_executable", type=Path)
    args = parser.parse_args()
    for path in (THREE_BODY, BIGSTICK):
        if not path.is_file():
            raise SystemExit(f"missing {path}")
    exe = str(args.imsrg_executable)
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        modelspace = pyIMSRG.ModelSpace(2, "O16", "O16")
        modelspace.SetHbarOmega(20.0)
        nn = root / "nn_O16.jcoupled64"
        pyIMSRG.ReadWrite().Write_jcoupled64(
            str(nn), vacuum_operator(modelspace, load_reference(REPOSITORY / "prototype/mrimsrg/data/O16_Nrefmax0_final")))
        common = [f"2bme={nn}", "fmt2=jcoupled64", f"3bme={THREE_BODY}", "3bme_type=no2b", "e3max=6",
                  "file3e1max=6", "file3e2max=12", "file3e3max=6", "reference=O16", "valence_space=O16", "A=16",
                  "hw=20", "emax=2", "BetaCM=0", "nucleon_mass_correction=false",
                  "denominator_partitioning=Epstein_Nesbet", "nsteps=1"]
        (root / "scratch").mkdir()
        log = run([exe, *common, "basis=HF", "method=magnus", "smax=0", "write_omega=true",
                   f"scratch={root / 'scratch'}", f"intfile={root / 'hf'}",
                   f"write_H_no2bpack={root / 'H_HF.no2bpack'}"], root / "hf.log")
        ehf = float(re.findall(r"^EHF = (\S+)", log, re.M)[-1])

        (root / "hf_basis.bin").symlink_to(root / "H_HF.no2bpack")
        bigstick(root, "n\nref\nauto\n2\n8 8\n0\n+\ny\n0\nno2\nhf_basis.bin\n20 0\nend\nex\n3\n", "ref.log")
        bigstick(root, "ba\nauto\n2\n8 8\n0\n+\ny\n0\nref\n", "ba.log")
        e_reference = first_energy(root / "ref.res")
        convert(root / "ref.bas", root / "ref.wfn", root / "H_HF.no2bpack", root / "ref.jref", nrefmax=0, hw=20.0,
                lambda3=False, basis_transformation=root / "hfC.mat", basis_orbits=root / "H_HF.no2bpack")

        mr = [exe, *common, "basis=oscillator", "method=magnus_adaptive", "core_generator=white-ncsm",
              f"mr_reference_file={root / 'ref.jref'}", "eta_criterion=1e-12", "mr_export_basis=natural"]
        log0 = run([*mr, "smax=0", f"flowfile={root / 'mr0.dat'}", f"intfile={root / 'mr0'}"], root / "mr0.log")
        e_mr0 = float(re.search(r"MR initial natural-basis Hamiltonian: E=(\S+)", log0).group(1))
        run([*mr, "smax=1", f"flowfile={root / 'mr1.dat'}", f"intfile={root / 'mr1'}",
             f"write_H_no2bpack={root / 'H_mr1.no2bpack'}"], root / "mr1.log")
        run([exe, *common, "basis=HF", "method=magnus_adaptive", "core_generator=white", "smax=1",
             "eta_criterion=1e-12", f"flowfile={root / 'sr1.dat'}", f"intfile={root / 'sr1'}"], root / "sr1.log")
        e_mr1, e_sr1 = last_flow_energy(root / "mr1.dat"), last_flow_energy(root / "sr1.dat")

        (root / "pd").mkdir()
        (root / "pd/h.bin").symlink_to(root / "H_mr1.no2bpack")
        bigstick(root / "pd", "n\npd\nauto\n2\n8 8\n0\n+\ny\n0\nno2\nh.bin\n20 0\nend\nex\n3\n", "pd.log")
        e_pd = first_energy(root / "pd/pd.res")

    print(f"HF {ehf:.7f}, BIGSTICK HF-basis reference {e_reference:.5f}, MR initial {e_mr0:.9f}")
    print(f"s=1: MR White-NCSM {e_mr1:.9f}, SR basis=HF White {e_sr1:.9f}, Nmax=0 of the MR-basis export {e_pd:.5f}")
    assert abs(e_mr0 - ehf) < 1e-6 and abs(e_reference - ehf) < 1e-4
    assert abs(e_mr1 - e_sr1) < 1e-8 and abs(e_pd - e_mr1) < 1e-4
    print("HF-basis reference chain is self-consistent")


if __name__ == "__main__":
    main()
