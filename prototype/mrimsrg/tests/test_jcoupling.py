from __future__ import annotations

import unittest

import numpy as np

from prototype.mrimsrg.jcoupling import (
    CoupledBlock,
    couple_scalar_two_body,
    coupled_scalar_contraction,
    extract_j_orbits,
    mr_lambda2_one_body_coupled,
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

    def test_mr_one_body_topologies_match_m_scheme(self) -> None:
        # Two radial copies for each species give nontrivial Hermitian and
        # anti-Hermitian matrices in every relevant pair channel.
        rows = []
        for index, (n, tz2) in enumerate(((0, -1), (0, +1), (1, -1), (1, +1))):
            for m2 in (-1, +1):
                rows.append([index, n, 0, 1, m2, tz2])
        orbits = np.asarray(rows, dtype=np.int32)
        layout = couple_scalar_two_body(np.zeros((8, 8, 8, 8)), orbits)
        rng = np.random.default_rng(314159)

        def random_blocks(hermiticity: int) -> tuple[CoupledBlock, ...]:
            result = []
            for block in layout:
                raw = rng.normal(size=block.matrix.shape)
                matrix = 0.5 * (raw + hermiticity * raw.T)
                result.append(
                    CoupledBlock(block.J, block.parity, block.Tz, block.pairs, matrix)
                )
            return tuple(result)

        x = random_blocks(-1)
        y = random_blocks(+1)
        lambda2 = random_blocks(+1)
        x_m = reconstruct_scalar_two_body(x, orbits)
        y_m = reconstruct_scalar_two_body(y, orbits)
        lambda_m = reconstruct_scalar_two_body(lambda2, orbits)

        direct = 0.25 * (
            np.einsum("iabc,deja,debc->ij", x_m, y_m, lambda_m, optimize=True)
            - np.einsum("iabc,deja,debc->ij", y_m, x_m, lambda_m, optimize=True)
        )
        crossed = np.einsum(
            "iabc,bejd,aecd->ij", x_m, y_m, lambda_m, optimize=True
        ) - np.einsum("iabc,bejd,aecd->ij", y_m, x_m, lambda_m, optimize=True)
        last = -0.5 * (
            np.einsum("iajb,cdae,cdbe->ij", x_m, y_m, lambda_m, optimize=True)
            - np.einsum("iajb,cdae,cdbe->ij", y_m, x_m, lambda_m, optimize=True)
        ) + 0.5 * (
            np.einsum("iajb,bcde,acde->ij", x_m, y_m, lambda_m, optimize=True)
            - np.einsum("iajb,bcde,acde->ij", y_m, x_m, lambda_m, optimize=True)
        )

        j_orbits = extract_j_orbits(orbits)

        def scalar_one_body_to_j(matrix: np.ndarray) -> np.ndarray:
            result = np.zeros((len(j_orbits), len(j_orbits)))
            for a, oa in enumerate(j_orbits):
                for b, ob in enumerate(j_orbits):
                    values = [
                        matrix[p, q]
                        for p in oa.substates
                        for q in ob.substates
                        if orbits[p, 4] == orbits[q, 4]
                    ]
                    if values:
                        result[a, b] = np.mean(values)
            return result

        coupled = mr_lambda2_one_body_coupled(x, y, lambda2, orbits, -1, +1)
        np.testing.assert_allclose(coupled.iv, scalar_one_body_to_j(direct), atol=2e-12, rtol=0)
        np.testing.assert_allclose(coupled.v, scalar_one_body_to_j(crossed), atol=2e-12, rtol=0)
        np.testing.assert_allclose(coupled.vi, scalar_one_body_to_j(last), atol=2e-12, rtol=0)


if __name__ == "__main__":
    unittest.main()
