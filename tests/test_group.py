from itertools import product

import numpy as np

from spgrep.group import (
    get_factor_system_from_little_group,
    get_little_group,
    is_matrix_group,
)
from spgrep.utils import NDArrayComplex, NDArrayInt, ndarray2d_to_integer_tuple


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


def check_cocycle_condition(
    rotations: NDArrayInt,
    factor_system: NDArrayComplex,
) -> bool:
    if not is_matrix_group(rotations):
        return False

    rotations_int = [ndarray2d_to_integer_tuple(r) for r in rotations]

    for i, j, k in product(range(len(rotations)), repeat=3):
        jk = rotations_int.index(ndarray2d_to_integer_tuple(rotations[j] @ rotations[k]))
        ij = rotations_int.index(ndarray2d_to_integer_tuple(rotations[i] @ rotations[j]))
        if not np.isclose(
            factor_system[i, jk] * factor_system[j, k], factor_system[ij, k] * factor_system[i, j]
        ):
            return False

    return True
