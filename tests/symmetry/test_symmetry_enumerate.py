from __future__ import annotations

import numpy as np
import pytest

from spgrep.rep.irreps import project_to_irrep
from spgrep.rep.representation import get_character, is_representation, is_unitary
from spgrep.symmetry.enumerate import enumerate_small_representations, enumerate_unitary_irreps
from spgrep.symmetry.group import (
    check_cocycle_condition,
    get_cayley_table,
    get_factor_system_from_little_group,
    get_little_group,
)
from spgrep.symmetry.representation import (
    check_spacegroup_representation,
    get_regular_representation,
)
from spgrep.utils import get_symmetry_from_hall_number


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


@pytest.mark.parametrize(
    "hall_number",
    [
        6,  # P2_1 (No. 4)
        350,  # P4_1 (No. 76)
        390,  # P-42_1m (No.113)
        431,  # P3_1 (No. 144)
        521,  # Pn-3m (No.224)
    ],
)
@pytest.mark.parametrize(
    "kpoint",
    [
        np.array([0.5, 0.5, 0]),
    ],
)
@pytest.mark.parametrize(
    "method",
    [
        "Neto",
        "random",
    ],
)
@pytest.mark.parametrize(
    "origin_shift",
    [
        np.array([1 / 3, 0, 0]),
        np.array([3 / 4, 0, 0]),
        np.array([3 / 4, 1 / 4, 0]),
    ],
)
def test_small_representation_with_origin_shift(hall_number, kpoint, method, origin_shift):
    rotations, translations = get_symmetry_from_hall_number(hall_number)
    little_rotations, little_translations, _ = get_little_group(rotations, translations, kpoint)
    assert len(little_rotations) > 0

    new_little_translations = []
    for R, tau in zip(little_rotations, little_translations):
        new_tau = np.remainder(tau + R @ origin_shift - origin_shift, 1)
        new_little_translations.append(new_tau)
    new_little_translations = np.array(new_little_translations)

    # Test factor system
    factor_system = get_factor_system_from_little_group(
        little_rotations, new_little_translations, kpoint
    )
    assert check_cocycle_condition(little_rotations, factor_system)

    # Test "weighted" irreps
    table = get_cayley_table(little_rotations)
    irreps, _ = enumerate_unitary_irreps(little_rotations, factor_system, method=method)
    assert sum(irrep.shape[1] ** 2 for irrep in irreps) == len(little_rotations)
    for irrep in irreps:
        assert is_representation(irrep, table, factor_system)

    # Test small representations
    small_reps, _ = enumerate_small_representations(
        little_rotations, new_little_translations, kpoint, method=method
    )
    assert sum(irrep.shape[1] ** 2 for irrep in small_reps) == len(little_rotations)
    for irrep in small_reps:
        assert check_spacegroup_representation(
            little_rotations, new_little_translations, kpoint, irrep
        )


def test_frobenius_schur_indicator(C4):
    irreps, indicators = enumerate_unitary_irreps(C4)
    assert sorted(indicators) == [0, 0, 1, 1]


def test_project_to_irrep(C3v):
    reg = get_regular_representation(C3v)
    irreps, _ = enumerate_unitary_irreps(C3v)

    count = 0
    for irrep in irreps:
        projected = project_to_irrep(reg, irrep)
        count += len(projected)

        for basis in projected:
            assert np.allclose(np.linalg.norm(basis, axis=1), 1)

    assert count == sum(irrep.shape[1] for irrep in irreps)
