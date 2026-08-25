from __future__ import annotations

import unittest

import numpy as np

from prototype.mrimsrg.jcoupling import (
    CoupledBlock,
    couple_scalar_two_body,
    coupled_scalar_contraction,
    reconstruct_scalar_two_body,
)


def _tiny_orbits() -> np.ndarray:
    # Two distinguishable j=1/2 multiplets.  This exercises canonical-pair
    # phases as well as identical-pair normalization without a large fixture.
    return np.asarray(
        [
            [0, 0, 0, 1, -1, -1],
            [0, 0, 0, 1, +1, -1],
            [1, 0, 0, 1, -1, +1],
            [1, 0, 0, 1, +1, +1],
        ],
        dtype=np.int32,
    )


class JCouplingTest(unittest.TestCase):
    def test_random_normalized_blocks_round_trip(self) -> None:
        orbits = _tiny_orbits()
        zero = np.zeros((4, 4, 4, 4))
        layout = couple_scalar_two_body(zero, orbits)
        rng = np.random.default_rng(918273)
        blocks = tuple(
            CoupledBlock(
                block.J,
                block.parity,
                block.Tz,
                block.pairs,
                (lambda raw: 0.5 * (raw + raw.T))(
                    rng.normal(size=block.matrix.shape)
                ),
            )
            for block in layout
        )
        dense = reconstruct_scalar_two_body(blocks, orbits)
        recovered = couple_scalar_two_body(dense, orbits)
        for expected, actual in zip(blocks, recovered):
            np.testing.assert_allclose(actual.matrix, expected.matrix, atol=2e-14, rtol=0)

    def test_coupled_contraction_matches_m_scheme(self) -> None:
        orbits = _tiny_orbits()
        layout = couple_scalar_two_body(np.zeros((4, 4, 4, 4)), orbits)
        rng = np.random.default_rng(271828)

        def random_blocks() -> tuple[CoupledBlock, ...]:
            return tuple(
                CoupledBlock(
                    block.J,
                    block.parity,
                    block.Tz,
                    block.pairs,
                    (lambda raw: 0.5 * (raw + raw.T))(
                        rng.normal(size=block.matrix.shape)
                    ),
                )
                for block in layout
            )

        left = random_blocks()
        right = random_blocks()
        left_m = reconstruct_scalar_two_body(left, orbits)
        right_m = reconstruct_scalar_two_body(right, orbits)
        expected = 0.25 * np.einsum("pqrs,pqrs->", left_m, right_m, optimize=True)
        self.assertAlmostEqual(coupled_scalar_contraction(left, right), expected, places=13)


if __name__ == "__main__":
    unittest.main()
