"""Representation-matrix related implementations for space group."""

from __future__ import annotations

from itertools import product

import numpy as np

from spgrep._constants import RTOL
from spgrep.symmetry.group import get_cayley_table
from spgrep.utils import (
    NDArrayComplex,
    NDArrayFloat,
    NDArrayInt,
    ndarray2d_to_integer_tuple,
)


def get_regular_representation(rotations: NDArrayInt) -> NDArrayInt:
    """Calculate regular representation of point group.

    Parameters
    ----------
    rotations: array, (order, 3, 3)

    Returns
    -------
    reg: array, (order, order, order)
        ``reg[k]`` is a representation matrix for ``rotations[k]``.
        If and only if ``np.dot(rotations[k], rotations[j]) == rotations[i]``, ``reg[k, i, j] == 1``.
    """
    n = len(rotations)
    table = get_cayley_table(rotations)

    reg = np.zeros((n, n, n), dtype=int)
    for k, j in product(range(n), repeat=2):
        reg[k, table[k, j], j] = 1

    return reg


def get_projective_regular_representation(
    rotations: NDArrayInt, factor_system: NDArrayComplex
) -> NDArrayComplex:
    """Calculate regular representation of space group with factor system.

    Parameters
    ----------
    rotations: array, (order, 3, 3)
    factor_system: array, (order, order)

    Returns
    -------
    reg: array, (order, order, order)
        ``reg[k]`` is a representation matrix for ``rotations[k]``.
        If and only if ``np.dot(rotations[k], rotations[j]) == rotations[i]``, ``reg[k, i, j] == factor_system[k, j]``.
    """
    n = len(rotations)
    table = get_cayley_table(rotations)

    reg = np.zeros((n, n, n), dtype=np.complex128)
    for k, j in product(range(n), repeat=2):
        reg[k, table[k, j], j] = factor_system[k, j]

    return reg


def check_spacegroup_representation(
    little_rotations: NDArrayInt,
    little_translations: NDArrayFloat,
    kpoint: NDArrayFloat,
    rep: NDArrayComplex,
    spinor_factor_system: NDArrayComplex | None = None,
    rtol: float = RTOL,
) -> bool:
    """Check definition of representation. This function works for primitive and conventional cell."""
    order = len(little_rotations)
    if spinor_factor_system is None:
        spinor_factor_system = np.ones((order, order), dtype=np.complex128)

    little_rotations_int = [ndarray2d_to_integer_tuple(rotation) for rotation in little_rotations]

    # Check if ``rep`` preserves multiplication
    for idx1, (r1, t1, m1) in enumerate(zip(little_rotations, little_translations, rep)):
        for idx2, (r2, t2, m2) in enumerate(zip(little_rotations, little_translations, rep)):
            r12 = r1 @ r2
            t12 = r1 @ t2 + t1
            idx = little_rotations_int.index(ndarray2d_to_integer_tuple(r12))
            # little_translations[idx] may differ from t12 by lattice translation.
            m12 = (
                spinor_factor_system[idx1, idx2]
                * np.exp(-2j * np.pi * np.dot(kpoint, t12 - little_translations[idx]))
                * rep[idx]
            )

            if not np.allclose(m12, m1 @ m2, rtol=rtol):
                return False

    return True
