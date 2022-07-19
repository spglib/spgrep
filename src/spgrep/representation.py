from __future__ import annotations

from itertools import product
from warnings import warn

import numpy as np

from spgrep.group import get_cayley_table
from spgrep.utils import NDArrayComplex, NDArrayInt


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

    reg = np.zeros((n, n, n), dtype=np.complex_)
    for k, j in product(range(n), repeat=2):
        reg[k, table[k, j], j] = factor_system[k, j]

    return reg


def get_intertwiner(
    rep1: NDArrayComplex,
    rep2: NDArrayComplex,
    rtol: float = 1e-5,
    max_num_random_generations: int = 4,
):
    """Calculate intertwiner matrix between ``rep1`` and ``rep2`` such that ``rep1 @ matrix == matrix @ rep2`` if they are equivalent.
    This function takes O(order * dim^4).

    Parameters
    ----------
    rep1: array, (order, dim, dim)
        Unitary irrep
    rep2: array, (order, dim, dim)
        Unitary irrep
    rtol: float
        Relative tolerance to distinguish difference eigenvalues
    max_num_random_generations: int
        Maximal number of trials to generate random matrix

    Returns
    -------
    matrix: array, (dim, dim)
    """
    assert rep1.shape == rep2.shape
    dim = rep1.shape[1]

    rng = np.random.default_rng(0)
    for _ in range(max_num_random_generations):
        random = rng.random((dim, dim)) + rng.random((dim, dim)) * 1j
        matrix = np.einsum("kil,lm,kjm->ij", rep1, random, np.conj(rep2))
        if not np.allclose(matrix, 0, rtol=rtol):
            return matrix

    warn("Failed to search all irreps. Try increasing max_num_random_generations.")
    return np.zeros((dim, dim))


def get_character(representation: NDArrayComplex) -> NDArrayComplex:
    """Calculate character of representation

    Parameters
    ----------
    representation: array, (order, dim, dim)

    Returns
    -------
    character: array, (order, )
    """
    character = np.einsum("ijj->i", representation, optimize="greedy")
    return character


def is_unitary(representation: NDArrayComplex) -> bool:
    dim = representation.shape[1]
    for matrix in representation:
        if not np.allclose(matrix @ np.conj(matrix.T), np.eye(dim)):
            return False
    return True


def is_projective_representation(
    rep: NDArrayComplex,
    table: NDArrayInt,
    factor_system: NDArrayComplex,
    rtol: float = 1e-5,
) -> bool:
    for i, ri in enumerate(rep):
        for j, rj in enumerate(rep):
            actual = ri @ rj
            expect = rep[table[i, j]] * factor_system[i, j]
            if not np.allclose(actual, expect, rtol=rtol):
                return False

    return True


def frobenius_schur_indicator(irrep: NDArrayComplex) -> int:
    """Inspect given unitary (projective) irrep is real, pseudo-real, or not unitary equivalent.

    .. math::
       \\mathrm{indicator} =
       \\frac{1}{|G|} \\sum_{ g \\in G } \\chi(g^{2})

    Parameters
    ----------
    irrep: array, (order, dim, dim)

    Returns
    -------
    indicator: int
        If indicator==1, it is real Reps.
        If indicator==-1, it is psedu-real Reps.
        Otherwise, it and adjoint Reps. are not equivalent.
    """
    order = irrep.shape[0]
    indicator = np.einsum("kij,kji->", irrep, irrep) / order
    indicator = int(np.around(np.real(indicator)))

    if indicator > 1:
        raise ValueError(f"Given representation is not irreducible: indicator={indicator}")

    return indicator
