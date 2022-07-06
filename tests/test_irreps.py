import numpy as np
import pytest

from spgrep.core import (
    get_spacegroup_irreps,
    get_spacegroup_irreps_from_primitive_symmetry,
)
from spgrep.irreps import get_character, get_irreps, get_regular_representation


def test_get_character(C3v):
    reg = get_regular_representation(C3v)
    actual = get_character(reg)
    expect = np.array([6, 0, 0, 0, 0, 0])
    assert np.allclose(actual, expect)


def test_get_irreps_C3v(C3v):
    reg = get_regular_representation(C3v)
    irreps = get_irreps(reg)
    # Check dimensions
    assert [irrep.shape[1] for irrep in irreps] == [1, 1, 2]
    # Check characters
    characters_expect = np.array(
        [
            [1, 1, 1, 1, 1, 1],  # A1
            [1, 1, 1, -1, -1, -1],  # A2
            [2, -1, -1, 0, 0, 0],  # E
        ]
    )
    characters_actual = np.array([get_character(irrep) for irrep in irreps])
    assert np.allclose(characters_actual, characters_expect)


@pytest.mark.skip
def test_get_spacegroup_irreps_from_primitive_symmetry(P42mnm):
    rotations, translations = P42mnm
    kpoint = np.array([0, 1 / 2, 0])  # X point
    get_spacegroup_irreps_from_primitive_symmetry(rotations, translations, kpoint)


@pytest.mark.skip
def test_get_irreps_Ia3d_H(Ia3d_H):
    pass


@pytest.mark.skip
def test_get_spacegroup_irreps(corundum_cell):
    lattice, positions, numbers = corundum_cell
    # kpoint: T for hR
    # kpoint = np.array([1 / 2, 1 / 2, -1 / 2])
    kpoint = np.array([0, 1, 1 / 2])
    irreps = get_spacegroup_irreps(*corundum_cell, kpoint=kpoint)  # noqa: F841
