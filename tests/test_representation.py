import numpy as np

from spgrep.irreps import enumerate_unitary_irreps, is_equivalent_irrep
from spgrep.representation import (
    get_character,
    get_intertwiner,
    get_regular_representation,
    project_to_irrep,
)


def test_get_character(C3v):
    reg = get_regular_representation(C3v)
    actual = get_character(reg)
    expect = np.array([6, 0, 0, 0, 0, 0])
    assert np.allclose(actual, expect)


def test_intertwiner():
    rep1 = np.array(
        [
            [[1.0 - 0.0j, 0.0 + 0.0j], [0.0 + 0.0j, 1.0 - 0.0j]],
            [[0.0 + 0.0j, 1.0 - 0.0j], [1.0 - 0.0j, 0.0 + 0.0j]],
            [[-0.0 - 1.0j, 0.0 + 0.0j], [0.0 + 0.0j, -0.0 + 1.0j]],
            [[0.0 + 0.0j, -0.0 - 1.0j], [-0.0 + 1.0j, 0.0 + 0.0j]],
        ]
    )
    rep2 = np.array(
        [
            [[1.0 - 0.0j, 0.0 + 0.0j], [0.0 + 0.0j, 1.0 - 0.0j]],
            [[0.0 + 0.0j, 1.0 - 0.0j], [1.0 - 0.0j, 0.0 + 0.0j]],
            [[-0.0 + 1.0j, 0.0 - 0.0j], [0.0 - 0.0j, 0.0 - 1.0j]],
            [[0.0 - 0.0j, -0.0 + 1.0j], [0.0 - 1.0j, 0.0 - 0.0j]],
        ]
    )

    intertwiner = get_intertwiner(rep1, rep2)
    assert is_equivalent_irrep(get_character(rep1), get_character(rep2))
    assert np.allclose(
        np.einsum("kil,lj->kij", rep1, intertwiner),
        np.einsum("il,klj->kij", intertwiner, rep2),
    )


def test_project_to_irrep(C3v):
    reg = get_regular_representation(C3v)
    irreps, _ = enumerate_unitary_irreps(C3v)

    count = 0
    for irrep in irreps:
        projected = project_to_irrep(reg, irrep)
        count += len(projected)

    assert count == sum(irrep.shape[1] for irrep in irreps)


def test_frobenius_schur_indicator(C4):
    irreps, indicators = enumerate_unitary_irreps(C4)
    assert sorted(indicators) == [0, 0, 1, 1]
