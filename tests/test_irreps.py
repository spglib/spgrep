import numpy as np
import pytest

from spgrep.core import (
    get_crystallographic_pointgroup_irreps_from_symmetry,
    get_spacegroup_irreps,
    get_spacegroup_irreps_from_primitive_symmetry,
)
from spgrep.group import get_cayley_table
from spgrep.irreps import enumerate_unitary_irreps_from_regular_representation
from spgrep.representation import (
    get_character,
    get_regular_representation,
    is_projective_representation,
    is_unitary,
)
from spgrep.transform import transform_symmetry_and_kpoint, unique_primitive_symmetry
from spgrep.utils import (
    NDArrayComplex,
    NDArrayFloat,
    NDArrayInt,
    ndarray2d_to_integer_tuple,
)


def test_get_irreps_random_C3v(C3v):
    reg = get_regular_representation(C3v)
    irreps = enumerate_unitary_irreps_from_regular_representation(
        reg=reg.astype(np.cdouble),
    )
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
def test_get_crystallographic_pointgroup_irreps_from_symmetry_C3v(method, C3v):
    irreps = get_crystallographic_pointgroup_irreps_from_symmetry(C3v, method=method)
    table = get_cayley_table(C3v)
    factor_system = np.ones((6, 6), dtype=np.complex_)
    assert np.sum([irrep.shape[1] ** 2 for irrep in irreps]) == 6
    for irrep in irreps:
        assert is_projective_representation(irrep, table, factor_system)


@pytest.mark.parametrize("method", [("Neto"), ("random")])
def test_get_spacegroup_irreps_from_primitive_symmetry_P42mnm(method, P42mnm):
    rotations, translations = P42mnm
    kpoint = np.array([0, 1 / 2, 0])  # X point
    irreps, mapping_little_group = get_spacegroup_irreps_from_primitive_symmetry(
        rotations, translations, kpoint, method=method
    )
    assert len(irreps) == 2
    assert [irrep.shape[1] for irrep in irreps] == [2, 2]

    little_rotations = rotations[mapping_little_group]
    little_translations = translations[mapping_little_group]
    for irrep in irreps:
        assert check_spacegroup_representation(
            little_rotations, little_translations, kpoint, irrep
        )
        assert is_unitary(irrep)


@pytest.mark.parametrize("method", [("Neto"), ("random")])
def test_get_spacegroup_irreps_from_primitive_symmetry_Ia3d(method, Ia3d):
    rotations, translations = Ia3d
    kpoint_conv = np.array([0, 1, 0])  # H point in conventional dual

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
    primitive_rotations, primitive_translations, mapping = unique_primitive_symmetry(
        primitive_rotations, primitive_translations
    )
    assert primitive_rotations.shape == (48, 3, 3)
    assert primitive_translations.shape == (48, 3)
    kpoint_prim = np.array([1 / 2, -1 / 2, 1 / 2])
    assert np.allclose(primitive_kpoint, kpoint_prim)

    primitive_irreps, mapping_little_group = get_spacegroup_irreps_from_primitive_symmetry(
        rotations=primitive_rotations,
        translations=primitive_translations,
        kpoint=primitive_kpoint,
        method=method,
    )
    # 48 = 2^2 + 2^2 + 2^2 + 6^2
    assert len(primitive_irreps) == 4
    assert sorted(irrep.shape[1] for irrep in primitive_irreps) == [2, 2, 2, 6]

    little_primitive_rotations = primitive_rotations[mapping_little_group]
    little_primitive_translations = primitive_translations[mapping_little_group]
    for irrep in primitive_irreps:
        assert check_spacegroup_representation(
            little_primitive_rotations, little_primitive_translations, primitive_kpoint, irrep
        )
        assert is_unitary(irrep)


@pytest.mark.parametrize("method", [("Neto"), ("random")])
def test_get_spacegroup_irreps(method, corundum_cell):
    # Corundum structure, R-3c (No. 167)
    # T point for hR
    kpoint = np.array([0, 1, 1 / 2])
    irreps, rotations, translations, mapping = get_spacegroup_irreps(
        *corundum_cell, kpoint=kpoint, method=method
    )
    # order=12, 12 = 2^2 + 2^2 + 2^2
    assert len(irreps) == 3
    assert [irrep.shape[1] for irrep in irreps] == [2, 2, 2]
    assert set(mapping) == set(range(36))  # order of little co-group in conventional

    little_rotations = rotations[mapping]
    little_translations = translations[mapping]
    for irrep in irreps:
        assert check_spacegroup_representation(
            little_rotations, little_translations, kpoint, irrep
        )
        assert is_unitary(irrep)


def check_spacegroup_representation(
    little_rotations: NDArrayInt,
    little_translations: NDArrayFloat,
    kpoint: NDArrayFloat,
    rep: NDArrayComplex,
):
    """Check definition of representation. This function works for primitive and conventional cell."""
    little_rotations_int = [ndarray2d_to_integer_tuple(rotation) for rotation in little_rotations]

    # Check if ``rep`` preserves multiplication
    for r1, t1, m1 in zip(little_rotations, little_translations, rep):
        for r2, t2, m2 in zip(little_rotations, little_translations, rep):
            r12 = r1 @ r2
            t12 = r1 @ t2 + t1
            idx = little_rotations_int.index(ndarray2d_to_integer_tuple(r12))
            # little_translations[idx] may differ from t12 by lattice translation.
            m12 = rep[idx] * np.exp(-2j * np.pi * np.dot(kpoint, t12 - little_translations[idx]))

            if not np.allclose(m12, m1 @ m2):
                return False

    return True
