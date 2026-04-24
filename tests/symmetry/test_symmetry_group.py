from __future__ import annotations

import numpy as np

from spgrep.symmetry.group import (
    check_cocycle_condition,
    decompose_by_maximal_space_subgroup,
    get_factor_system_from_little_group,
    get_little_group,
    get_little_group_of_pm_k,
    is_matrix_group,
)
from spgrep.utils import get_symmetry_from_hall_number


def test_is_matrix_group(C3v):
    assert is_matrix_group(C3v)


def test_get_little_group_and_factor_system(P42mnm):
    rotations, translations = P42mnm
    kpoint = np.array([0, 1 / 2, 0])  # X point

    # Check little group
    little_rotations, little_translations, mapping_little_group = get_little_group(
        rotations, translations, kpoint
    )
    assert little_rotations.shape == (8, 3, 3)
    assert little_translations.shape == (8, 3)
    assert len(mapping_little_group) == 8

    # Check factor system
    factor_system = get_factor_system_from_little_group(
        little_rotations, little_translations, kpoint
    )
    assert check_cocycle_condition(little_rotations, factor_system)


def test_get_little_group_of_pm_k_coincides_with_little_group_when_2k_is_zero(P42mnm):
    # X = (0, 1/2, 0): 2k is a reciprocal lattice vector, so G^{k,k-bar} = G^k.
    rotations, translations = P42mnm
    kpoint = np.array([0, 1 / 2, 0])

    little_rotations, little_translations, mapping_little_group = get_little_group(
        rotations, translations, kpoint
    )
    pm_rotations, pm_translations, pm_mapping, flip_k = get_little_group_of_pm_k(
        rotations, translations, kpoint
    )

    np.testing.assert_array_equal(pm_rotations, little_rotations)
    np.testing.assert_allclose(pm_translations, little_translations)
    np.testing.assert_array_equal(pm_mapping, mapping_little_group)
    assert flip_k.shape == (len(little_rotations),)
    assert flip_k.dtype == np.bool_
    assert not np.any(flip_k)


def test_get_little_group_of_pm_k_on_pm3m_lambda():
    # SG 221 (Pm-3m), Hall 517. Along Lambda k = (u, u, u) we have
    # 2k = (2u, 2u, 2u) != 0 mod 1 for generic u, so G^{k,k-bar} strictly
    # extends G^k. Little cogroup of k is D3 along [111] (order 6); little
    # cogroup of {k, -k} is D3d (order 12) with the extra 6 ops mapping k
    # to -k (see the Bilbao Lambda-line listing for Pm-3m).
    rotations, translations = get_symmetry_from_hall_number(hall_number=517)
    kpoint = np.array([1 / 4, 1 / 4, 1 / 4])

    pm_rotations, pm_translations, pm_mapping, flip_k = get_little_group_of_pm_k(
        rotations, translations, kpoint
    )
    assert pm_rotations.shape == (12, 3, 3)
    assert pm_translations.shape == (12, 3)
    assert pm_mapping.shape == (12,)
    assert flip_k.shape == (12,)
    assert flip_k.dtype == np.bool_
    assert np.count_nonzero(~flip_k) == 6
    assert np.count_nonzero(flip_k) == 6

    # Stabilizing ops match get_little_group.
    little_rotations, _, _ = get_little_group(rotations, translations, kpoint)
    np.testing.assert_array_equal(pm_rotations[~flip_k], little_rotations)

    # k-flipping ops satisfy S^T k == -k (mod 1).
    for rot in pm_rotations[flip_k]:
        residual = rot.T @ kpoint + kpoint
        residual -= np.rint(residual)
        np.testing.assert_allclose(residual, 0, atol=1e-8)


def test_decompose_by_maximal_space_subgroup(P42mnm_type3):
    rotations, translations, time_reversals, _ = P42mnm_type3
    (xsg_indices, time_reversal_indices, _) = decompose_by_maximal_space_subgroup(
        rotations, translations, time_reversals
    )
    assert len(xsg_indices) == 8
    assert len(time_reversal_indices) == 8
    assert sorted(xsg_indices + time_reversal_indices) == list(range(16))
