"""Representation-matrix related implementations."""

from __future__ import annotations

from warnings import warn

import numpy as np

from spgrep._constants import ATOL, MAX_NUM_RANDOM_GENERATIONS, RTOL
from spgrep.utils import (
    NDArrayComplex,
    NDArrayFloat,
    NDArrayInt,
)


def get_intertwiner(
    rep1: NDArrayComplex,
    rep2: NDArrayComplex,
    atol: float = ATOL,
    max_num_random_generations: int = MAX_NUM_RANDOM_GENERATIONS,
) -> NDArrayComplex:
    """Calculate intertwiner matrix between ``rep1`` and ``rep2`` such that ``rep1 @ matrix == matrix @ rep2`` if they are equivalent.

    The determinant of ``matrix`` is scaled to be unity.

    This function takes O(order * dim^4).

    Parameters
    ----------
    rep1: array, (order, dim, dim)
        Unitary irrep
    rep2: array, (order, dim, dim)
        Unitary irrep
    atol: float
        Absolute tolerance to distinguish difference eigenvalues
    max_num_random_generations: int
        Maximum number of trials to generate random matrix

    Returns
    -------
    matrix: array, (dim, dim)
    """
    assert rep1.shape == rep2.shape
    dim = rep1.shape[1]

    rng = np.random.default_rng(0)
    for _ in range(max_num_random_generations):
        random = rng.random((dim, dim)) + rng.random((dim, dim)) * 1j
        matrix: NDArrayComplex = np.einsum(
            "kil,lm,kjm->ij", rep1, random, np.conj(rep2), optimize="greedy"
        )
        if not np.allclose(matrix, 0, atol=atol):
            # Scale such that determinant is unity
            matrix /= np.linalg.det(matrix) ** (1 / dim)
            return matrix

    warn("Failed to search all irreps. Try increasing max_num_random_generations.")
    return np.zeros((dim, dim), dtype=np.complex128)


def get_character(representation: NDArrayComplex) -> NDArrayComplex:
    """Calculate character of representation.

    Parameters
    ----------
    representation: array, (order, dim, dim)

    Returns
    -------
    character: array, (order, )
    """
    character = np.einsum("ijj->i", representation, optimize="greedy").astype(np.complex128)
    return character


def is_unitary(representation: NDArrayComplex) -> bool:
    """Return true if given representation is unitary."""
    dim = representation.shape[1]
    for matrix in representation:
        if not np.allclose(matrix @ np.conj(matrix.T), np.eye(dim)):
            return False
    return True


def is_representation(
    rep: NDArrayComplex,
    table: NDArrayInt,
    factor_system: NDArrayComplex | None = None,
    rtol: float = RTOL,
    atol: float = ATOL,
) -> bool:
    """Return true if given matrix function is a (projective) representation with given factor system."""
    order = rep.shape[0]
    if factor_system is None:
        factor_system = np.ones((order, order), dtype=np.complex128)

    for i, ri in enumerate(rep):
        for j, rj in enumerate(rep):
            actual = ri @ rj
            expect = rep[table[i, j]] * factor_system[i, j]
            if not np.allclose(actual, expect, rtol=rtol, atol=atol):
                return False

    return True


def get_direct_product(
    rep1: NDArrayComplex | NDArrayFloat, rep2: NDArrayComplex | NDArrayFloat
) -> NDArrayComplex | NDArrayFloat:
    """Return Knocker product of two representations.

    Parameters
    ----------
    rep1: array, (order, dim1, dim1)
    rep2: array, (order, dim2, dim2)

    Returns
    -------
    direct: (order, dim1 * dim2, dim1 * dim2)
    """
    order = rep1.shape[0]
    dim1 = rep1.shape[1]
    dim2 = rep2.shape[1]

    if rep1.shape != (order, dim1, dim1) or rep2.shape != (order, dim2, dim2):
        raise ValueError("Inconsistent shapes.")

    direct = (rep1[:, :, None, :, None] * rep2[:, None, :, None, :]).reshape(
        order, dim1 * dim2, dim1 * dim2
    )
    return direct
