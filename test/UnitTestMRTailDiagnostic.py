#!/usr/bin/env python3

"""Exercise the White-NCSM channel diagnostic on the real He4 fixture."""

from pathlib import Path
import sys
import tempfile

import pyIMSRG


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))
sys.path.insert(0, str(REPOSITORY / "test"))

from UnitTestMRDriver import vacuum_operator  # noqa: E402
from prototype.mrimsrg.diagnose_white_ncsm_tail import diagnose_stage  # noqa: E402
from prototype.mrimsrg.export_jref import export_reference  # noqa: E402
from prototype.mrimsrg.reference_io import load_reference  # noqa: E402


def main() -> None:
    fixture = REPOSITORY / "prototype/mrimsrg/data/He4_Nrefmax2"
    data = load_reference(fixture)
    with tempfile.TemporaryDirectory() as temporary_directory:
        root = Path(temporary_directory)
        reference_path = root / "He4_Nrefmax2.jref"
        hamiltonian_path = root / "He4_bare.jcoupled64"
        export_reference(fixture, reference_path)

        vacuum_modelspace = pyIMSRG.ModelSpace(2, "He4", "He4")
        vacuum_modelspace.SetHbarOmega(20.0)
        pyIMSRG.ReadWrite().Write_jcoupled64(
            str(hamiltonian_path), vacuum_operator(vacuum_modelspace, data)
        )

        modelspace = pyIMSRG.ModelSpace(2, "He4", "He4")
        modelspace.SetHbarOmega(20.0)
        modelspace.SetReferenceOcc(
            pyIMSRG.MRReference.ReadOccupationMap(
                modelspace, str(reference_path)
            )
        )
        reference = pyIMSRG.MRReference.ReadBinary(
            modelspace, str(reference_path)
        )
        stage = diagnose_stage(
            pyIMSRG,
            modelspace,
            reference,
            "fixture",
            hamiltonian_path,
            top=8,
        )

        assert stage["selected_channel_count"] > 0
        assert stage["nonzero_eta_channel_count"] > 0
        assert stage["eta_norm"] > 0.0
        assert abs(stage["eta_norm"] - stage["eta_norm_reconstructed"]) < 1e-12
        assert stage["white_ncsm_numerator_norm_mev"] > 0.0
        assert stage["selected_unweighted_hamiltonian_norm_mev"] > 0.0
        assert 0.0 < stage["top_channel_eta_norm_squared_fraction"] <= 1.0
        assert len(stage["top_channels"]) == 8
        contributions = [
            entry["eta_norm_squared_contribution"]
            for entry in stage["top_channels"]
        ]
        assert contributions == sorted(contributions, reverse=True)
        for entry in stage["top_channels"]:
            reconstructed = entry["hamiltonian_mev"] * (
                entry["forward_weight"] / entry["forward_denominator_mev"]
                - entry["reverse_weight"] / entry["reverse_denominator_mev"]
            )
            assert abs(reconstructed - entry["eta_mev"]) < 1e-13

    print("MR White-NCSM tail diagnostic test passed")


if __name__ == "__main__":
    main()
