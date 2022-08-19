import numpy as np

from spgrep.group import (
    check_cocycle_condition,
    get_factor_system_from_little_group,
    get_little_group,
    is_matrix_group,
)


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
