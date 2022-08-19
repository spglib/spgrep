"""Group-theory related functions."""
from __future__ import annotations

from itertools import product

import numpy as np

from spgrep.utils import (
    NDArrayComplex,
    NDArrayFloat,
    NDArrayInt,
    ndarray2d_to_integer_tuple,
)


def get_cayley_table(rotations: NDArrayInt) -> NDArrayInt:
    """Calculate Group multiplication table.

    Parameters
    ----------
    rotations: (order, 3, 3)

    Returns
    -------
    table: (order, order)
        ``table[i, j] = k`` if ``rotations[i] @ rotations[j] == rotations[k]``
    """
    rotations_list = [ndarray2d_to_integer_tuple(r) for r in rotations]

    order = rotations.shape[0]
    table = [[-1 for _ in range(order)] for _ in range(order)]
    for i, gi in enumerate(rotations):
        for j, gj in enumerate(rotations):
            gk = gi @ gj
            k = rotations_list.index(ndarray2d_to_integer_tuple(gk))
            if table[i][j] != -1:
                ValueError("Should specify a matrix group.")
            table[i][j] = k

    return np.array(table)


def get_identity_index(table: NDArrayInt) -> int:
    """Return index for identity of group."""
    order = table.shape[0]
    for i in range(order):
        if np.all(table[i, :] == np.arange(order)):
            return i

    raise ValueError("Unreachable!")


def get_inverse_index(table: NDArrayInt, idx: int) -> int:
    """Return index of inverse of ``idx`` element in ``table``."""
    order = table.shape[0]
    id_idx = get_identity_index(table)
    for i in range(order):
        if table[idx, i] == id_idx:
            return i

    raise ValueError("Unreachable!")


def get_order(table: NDArrayInt, idx: int) -> int:
    """Return order of element ``idx`` in ``table``."""
    id_idx = get_identity_index(table)
    ret = 1
    tmp = idx
    while tmp != id_idx:
        tmp = table[tmp, idx]
        ret += 1
    return ret


def is_matrix_group(rotations: NDArrayInt) -> bool:
    """Return True iff given integer matrices forms group."""
    try:
        table = get_cayley_table(rotations)
    except ValueError:
        return False

    # Check if each element appears only once in a row
    order = rotations.shape[0]
    for i in range(order):
        if set(table[i]) != set(range(order)):
            return False

    return True


def get_factor_system_from_little_group(
    little_rotations: NDArrayInt,
    little_translations: NDArrayFloat,
    kpoint: NDArrayFloat,
) -> NDArrayComplex:
    r"""Calculate factor system of projective representation of little co-group.

    .. math::
       D^{\mathbf{k}}_{p}(S_{i}) D^{\mathbf{k}}_{p}(S_{j})
       = \exp \left( -i \mathbf{g}_{i} \cdot \mathbf{w}_{j} \right) D^{\mathbf{k}}_{p}(S_{k})

    where :math:`S_{i}S_{j} = S_{k}` and :math:`\mathbf{g}_{i} = S_{i}^{-1} \mathbf{k} - \mathbf{k}`.

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
        Factor system of representations of little co-group that have one-to-one correspondence to small representations
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
    atol: float = 1e-8,
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
        ``(rotations[i], translations[i])`` belongs to the little group of given space space group and kpoint.
    """
    little_rotations = []
    little_translations = []
    mapping_little_group = []

    for i, (rotation, translation) in enumerate(zip(rotations, translations)):
        residual = rotation.T @ kpoint - kpoint
        residual = residual - np.rint(residual)
        if not np.allclose(residual, 0, atol=atol):
            continue
        little_rotations.append(rotation)
        little_translations.append(translation)
        mapping_little_group.append(i)

    return (
        np.array(little_rotations),
        np.array(little_translations),
        np.array(mapping_little_group),
    )


def check_cocycle_condition(
    rotations: NDArrayInt,
    factor_system: NDArrayComplex,
) -> bool:
    """Return true if given factor system satisfies the cocycle condition."""
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
