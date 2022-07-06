from __future__ import annotations

from itertools import product

import numpy as np

from spgrep.utils import (
    NDArrayComplex,
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


def get_factor_system_from_little_group(
    little_rotations: NDArrayInt,
    little_translations: NDArrayFloat,
    kpoint: NDArrayFloat,
) -> NDArrayComplex:
    """Calculate factor system of projective representation of given little group:

    .. math::
       D^{\\mathbf{k}}_{p}(S_{i}) D^{\\mathbf{k}}_{p}(S_{j})
       = \\exp \\left( -i \\mathbf{g}_{i} \\cdot \\mathbf{w}_{j} \\right) D^{\\mathbf{k}}_{p}(S_{k})

    where :math:`S_{i}S_{j} = S_{k}` and :math:`\\mathbf{g}_{i} = S_{i}^{-1} \\mathbf{k} - \\mathbf{k}`.

    Parameters
    ----------
    little_rotations: array, (order, 3, 3)
        Linear parts of coset of little group stabilizing ``kpoint``.
    little_translations: array, (order, 3)
        Translation parts of coset of little group stabilizing ``kpoint``.
    kpoint: array, (3, )

    Returns
    -------
    factor_system: array, (order, order)
        Factor system of small representation of given space group and kpoint.
    """
    n = len(little_rotations)

    residuals = np.zeros((n, 3))
    for i, rotation in enumerate(little_rotations):
        # Never take modulus!
        residuals[i] = rotation.T @ kpoint - kpoint

    factor_system = np.zeros((n, n), dtype=np.complex_)
    for i, residual in enumerate(residuals):
        for j, translation in enumerate(little_translations):
            factor_system[i, j] = np.exp(-2j * np.pi * np.dot(residual, translation))

    return factor_system


def get_little_group(
    rotations: NDArrayInt,
    translations: NDArrayFloat,
    kpoint: NDArrayFloat,
    rtol: float = 1e-5,
) -> tuple[NDArrayInt, NDArrayFloat, NDArrayInt]:
    """Return coset of little group of given space group which stabilize kpoint under rotations.

    Parameters
    ----------
    rotations: array, (order, 3, 3)
    translations: array, (order, 3)
    kpoint: array, (3, )

    Returns
    -------
    little_rotations: array, (little_group_order, 3, 3)
    little_translations: array, (little_group_order, 3)
    mapping_little_group: array, (little_group_order, )
        Let ``i = mapping_little_group[idx]``.
        (rotations[i], translations[i]) belongs to the little group of given space space group and kpoint.
    """
    little_rotations = []
    little_translations = []
    mapping_little_group = []

    for i, (rotation, translation) in enumerate(zip(rotations, translations)):
        residual = rotation.T @ kpoint - kpoint
        residual = residual - np.rint(residual)
        if not np.allclose(residual, 0, rtol=rtol):
            continue
        little_rotations.append(rotation)
        little_translations.append(translation)
        mapping_little_group.append(i)

    return (
        np.array(little_rotations),
        np.array(little_translations),
        np.array(mapping_little_group),
    )
