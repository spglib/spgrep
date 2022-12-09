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


def get_cayley_table(
    rotations: NDArrayInt, time_reversals: NDArrayInt | None = None
) -> NDArrayInt:
    """Calculate Group multiplication table.

    Parameters
    ----------
    rotations: array[int], (order, 3, 3)
    time_reversals: (Optional) array[int], (order, )

    Returns
    -------
    table: (order, order)
        ``table[i, j] = k`` if ``rotations[i] @ rotations[j] == rotations[k]``
    """
    order = rotations.shape[0]
    if time_reversals is None:
        time_reversals = np.zeros((order,), dtype=np.int_)

    operations_list = [
        (ndarray2d_to_integer_tuple(r), tr) for r, tr in zip(rotations, time_reversals)
    ]

    table = [[-1 for _ in range(order)] for _ in range(order)]
    for i, (ri, tri) in enumerate(zip(rotations, time_reversals)):
        for j, (rj, trj) in enumerate(zip(rotations, time_reversals)):
            rk = ri @ rj
            trk = tri != trj
            k = operations_list.index((ndarray2d_to_integer_tuple(rk), trk))
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


def decompose_by_maximal_space_subgroup(
    rotations: NDArrayInt,
    translations: NDArrayFloat,
    time_reversals: NDArrayInt,
) -> tuple[list[int], list[int], NDArrayInt, NDArrayFloat] | None:
    r"""Coset-decompose magnetic space group :math:`M` by its maximal space subgroup (XSG) :math:`D(M)`.

    If given magnetic space group is type I, return None.

    .. math::
        M = D(M) \sqcup (\mathbf{R}_{0}, \mathbf{v}_{0})1' D(M)

    Returns
    -------
    xsg_indices: list[int]
        List of indices for XSG
    time_reversal_indices: list[int]
        Let ``xsg_indices[i]`` = :math:`(\mathbf{W}_{i}, \mathbf{w}_{i})`.
        Then, ``time_reversal_indices[i]`` :math:`\equiv (\mathbf{R}_{0}, \mathbf{v}_{0}) (\mathbf{W}_{i}, \mathbf{w}_{i}) 1'`.
    conjugator_rotation: array[int], (3, 3)
        Rotation part of the coset representative :math:`\mathbf{R}_{0}`
    conjugator_translation: array, (3, )
        Translation part of the coset representative :math:`\mathbf{v}_{0}`
    """
    if np.all(time_reversals == 1):
        # Type-I MSG
        return None

    # Search for coset representative
    conjugator_rotation = None
    conjugator_translation = None
    for rot, trans, tr in zip(rotations, translations, time_reversals):
        if tr == 0:
            continue
        conjugator_rotation = rot.copy()
        conjugator_translation = trans.copy()
        break

    # Map XSG to time-reversal operations
    order = len(rotations)
    xsg_indices = list(np.arange(order)[time_reversals == 0])
    time_reversal_indices = []
    for i in xsg_indices:
        rot = rotations[i]
        trans = translations[i]
        # (R0, v0) (R, v) = (R0 @ R, R0 @ v + v0)
        new_rot = conjugator_rotation @ rot
        new_trans = conjugator_rotation @ trans + conjugator_translation
        for j in range(order):
            if time_reversals[j] == 0:
                continue
            if not np.allclose(new_rot, rotations[j]):
                continue
            diff = np.remainder(new_trans - translations[j], 1)
            diff -= np.rint(diff)
            if np.allclose(diff, 0):
                time_reversal_indices.append(j)
                break

    return xsg_indices, time_reversal_indices, conjugator_rotation, conjugator_translation
