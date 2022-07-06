from __future__ import annotations

from itertools import product

import numpy as np

from spgrep.utils import (
    NDArrayFloat,
    NDArrayInt,
    is_integer_array,
    ndarray2d_to_integer_tuple,
)


def is_matrix_group(rotations: NDArrayInt) -> bool:
    """
    Return True iff given integer matrices forms group.
    """
    rotations_set = {ndarray2d_to_integer_tuple(r) for r in rotations}

    for g1, g2 in product(rotations, repeat=2):
        # g1^-1 * g2 should be included in rotations_set
        g1_inv = np.linalg.inv(g1)
        if not is_integer_array(g1_inv):
            return False
        g1_inv = np.around(g1_inv).astype(int)
        g1_inv_g2 = np.dot(g1_inv, g2)
        if ndarray2d_to_integer_tuple(g1_inv_g2) not in rotations_set:
            return False

    return True


def get_factor_system(rotations: NDArrayInt, translations: NDArrayFloat, kpoint: NDArrayFloat):
    """


    Parameters
    ----------
    rotations: array, (order, 3, 3)
        Assume a fractional coordinates `x` are transformed by the i-th symmetry operation as follows:
            np.dot(rotations[i, :, :], x) + translations[i, :]
    translations: array, (order, 3)
    kpoint: array, (3, )

    Returns
    -------
    factor_system: array, (order, order)
        Factor system of small representation of given space group and kpoint.
    """
    raise NotImplementedError


def get_little_group(
    rotations: NDArrayInt,
    translations: NDArrayFloat,
    kpoint: NDArrayFloat,
    rtol: float = 1e-5,
) -> tuple[NDArrayInt, NDArrayFloat]:
    """Return coset of little group of given space group which stabilize kpoint under rotations."""
    little_rotations = []
    little_translations = []

    for rotation, translation in zip(rotations, translations):
        residual = rotation.T @ kpoint - kpoint
        residual = residual - np.rint(residual)
        if not np.allclose(residual, 0, rtol=rtol):
            continue
        little_rotations.append(rotation)
        little_translations.append(translation)

    return little_rotations, little_translations
