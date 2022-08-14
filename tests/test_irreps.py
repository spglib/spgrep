from __future__ import annotations

from itertools import product

import numpy as np
import pytest

from spgrep.core import (
    get_crystallographic_pointgroup_irreps_from_symmetry,
    get_spacegroup_irreps,
    get_spacegroup_irreps_from_primitive_symmetry,
)
from spgrep.group import get_cayley_table
from spgrep.irreps import enumerate_unitary_irreps, is_equivalent_irrep
from spgrep.pointgroup import pg_dataset
from spgrep.representation import (
    check_spacegroup_representation,
    get_character,
    is_representation,
    is_unitary,
)
from spgrep.transform import transform_symmetry_and_kpoint, unique_primitive_symmetry
from spgrep.utils import NDArrayComplex


@pytest.mark.parametrize("method", [("Neto"), ("random")])
def test_get_irreps_random_C3v(method, C3v):
    irreps, _ = enumerate_unitary_irreps(C3v, method=method)

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

    for irrep in irreps:
        assert is_unitary(irrep)


@pytest.mark.parametrize("method", [("Neto"), ("random")])
def test_get_crystallographic_pointgroup_irreps(method):
    for pg_symbol, groups in pg_dataset.items():
        for rotations in groups:
            if pg_symbol != "4":
                continue
            print(pg_symbol)
            irreps = get_crystallographic_pointgroup_irreps_from_symmetry(
                np.array(rotations), method=method
            )

            # Check exhaustiveness
            order = len(rotations)
            assert np.sum([irrep.shape[1] ** 2 for irrep in irreps]) == order

            # Check representation's property
            table = get_cayley_table(np.array(rotations))
            for irrep in irreps:
                assert is_representation(irrep, table)

            assert is_unique_irreps(irreps)


@pytest.mark.parametrize("method", [("Neto"), ("random")])
@pytest.mark.parametrize(
    "kpoint,shape_expect",
    [
        (np.array([0, 1 / 2, 0]), [2, 2]),  # X point
        (np.array([0, 0, 1 / 2]), [2, 2, 2, 2]),  # Z point
    ],
)
def test_get_spacegroup_irreps_from_primitive_symmetry_P42mnm(
    method, kpoint, shape_expect, P42mnm
):
    rotations, translations = P42mnm
    irreps, mapping_little_group = get_spacegroup_irreps_from_primitive_symmetry(
        rotations, translations, kpoint, method=method
    )
    assert [irrep.shape[1] for irrep in irreps] == shape_expect

    little_rotations = rotations[mapping_little_group]
    little_translations = translations[mapping_little_group]
    for irrep in irreps:
        assert check_spacegroup_representation(
            little_rotations, little_translations, kpoint, irrep
        )
        assert is_unitary(irrep)

    assert is_unique_irreps(irreps)


@pytest.mark.parametrize("method", [("Neto"), ("random")])
@pytest.mark.parametrize(
    "kpoint_conv,kpoint_prim_expect,shape_expect",
    [
        (np.array([0, 1, 0]), np.array([1 / 2, -1 / 2, 1 / 2]), [2, 2, 2, 6]),  # H point
        (np.array([1 / 2, 1 / 2, 0]), np.array([0, 0, 1 / 2]), [2, 2]),  # N point
    ],
)
def test_get_spacegroup_irreps_from_primitive_symmetry_Ia3d(
    method, kpoint_conv, kpoint_prim_expect, shape_expect, Ia3d
):
    rotations, translations = Ia3d

    # TODO: Refactor to function
    # Transform to primitive
    to_primitive = np.array(
        [
            [-1 / 2, 1 / 2, 1 / 2],
            [1 / 2, -1 / 2, 1 / 2],
            [1 / 2, 1 / 2, -1 / 2],
        ]
    )
    primitive_rotations, primitive_translations, primitive_kpoint = transform_symmetry_and_kpoint(
        to_primitive, rotations, translations, kpoint_conv
    )
    primitive_rotations, primitive_translations, _ = unique_primitive_symmetry(
        primitive_rotations, primitive_translations
    )
    assert primitive_rotations.shape == (48, 3, 3)
    assert primitive_translations.shape == (48, 3)
    assert np.allclose(primitive_kpoint, kpoint_prim_expect)

    primitive_irreps, mapping_little_group = get_spacegroup_irreps_from_primitive_symmetry(
        rotations=primitive_rotations,
        translations=primitive_translations,
        kpoint=primitive_kpoint,
        method=method,
    )
    assert np.sum([irrep.shape[1] ** 2 for irrep in primitive_irreps]) == len(mapping_little_group)
    assert sorted(irrep.shape[1] for irrep in primitive_irreps) == shape_expect

    little_primitive_rotations = primitive_rotations[mapping_little_group]
    little_primitive_translations = primitive_translations[mapping_little_group]
    for irrep in primitive_irreps:
        assert check_spacegroup_representation(
            little_primitive_rotations, little_primitive_translations, primitive_kpoint, irrep
        )
        assert is_unitary(irrep)

    assert is_unique_irreps(primitive_irreps)


@pytest.mark.parametrize("method", [("Neto"), ("random")])
@pytest.mark.parametrize(
    "kpoint,shape_expect,num_sym_expect",
    [
        (np.array([0, 1, 1 / 2]), [2, 2, 2], 36),  # T point for hR
        (np.array([-1 / 2, 1 / 2, 1 / 2]), [2], 12),  # L point for hR
    ],
)
def test_get_spacegroup_irreps(method, kpoint, shape_expect, num_sym_expect, corundum_cell):
    # Corundum structure, R-3c (No. 167)
    irreps, rotations, translations, mapping = get_spacegroup_irreps(
        *corundum_cell, kpoint=kpoint, method=method
    )
    assert [irrep.shape[1] for irrep in irreps] == shape_expect
    assert len(mapping) == num_sym_expect  # order of little co-group in conventional

    little_rotations = rotations[mapping]
    little_translations = translations[mapping]
    for irrep in irreps:
        assert check_spacegroup_representation(
            little_rotations, little_translations, kpoint, irrep
        )
        assert is_unitary(irrep)

    assert is_unique_irreps(irreps)


def is_unique_irreps(irreps: list[NDArrayComplex]):
    characters = [get_character(irrep) for irrep in irreps]
    for (i, ci), (j, cj) in product(enumerate(characters), repeat=2):
        if is_equivalent_irrep(ci, cj) != (i == j):
            return False
    return True


def test_tetragonal():
    # -42m
    rotations = np.array(
        [
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            [[0, -1, 0], [-1, 0, 0], [1, 1, 1]],
            [[-1, -1, -1], [0, 0, 1], [0, 1, 0]],
            [[0, 0, -1], [1, 1, 1], [-1, 0, 0]],
            [[0, 0, 1], [1, 0, 0], [-1, -1, -1]],
            [[1, 1, 1], [0, -1, 0], [0, 0, -1]],
            [[-1, 0, 0], [0, 0, -1], [0, -1, 0]],
            [[0, 1, 0], [-1, -1, -1], [1, 0, 0]],
        ]
    )
    irreps, _ = enumerate_unitary_irreps(rotations, method="Neto")
    assert len(irreps) == 5
