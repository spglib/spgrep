import numpy as np

from spgrep.group import get_little_group, is_matrix_group


def test_is_matrix_group(C3v):
    assert is_matrix_group(C3v)


def test_get_little_group(P3m1):
    rotations, translations = P3m1
    kpoint = np.array([1 / 2, 0, 0])  # M point

    little_rotations, little_translations = get_little_group(rotations, translations, kpoint)
    little_rotations_expect = np.array(
        [
            [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
            ],
            [
                [1, 0, 0],
                [1, -1, 0],
                [0, 0, 1],
            ],
        ]
    )
    little_translations_expect = np.array(
        [
            [0, 0, 0],
            [0, 0, 0],
        ]
    )
    assert np.allclose(little_rotations, little_rotations_expect)
    assert np.allclose(little_translations, little_translations_expect)
