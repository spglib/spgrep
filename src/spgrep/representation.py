from __future__ import annotations

from itertools import product
from warnings import warn

import numpy as np

from spgrep.group import get_cayley_table
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

    reg = np.zeros((n, n, n), dtype=np.complex_)
    for k, j in product(range(n), repeat=2):
        reg[k, table[k, j], j] = factor_system[k, j]

    return reg


def get_intertwiner(
    rep1: NDArrayComplex,
    rep2: NDArrayComplex,
    atol: float = 1e-8,
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
    atol: float
        Absolute tolerance to distinguish difference eigenvalues
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
        matrix = np.einsum("kil,lm,kjm->ij", rep1, random, np.conj(rep2), optimize="greedy")
        if not np.allclose(matrix, 0, atol=atol):
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
    character = np.einsum("ijj->i", representation, optimize="greedy").astype(np.complex_)
    return character


def project_to_irrep(
    representation: NDArrayComplex,
    irrep: NDArrayComplex,
    atol: float = 1e-8,
) -> list[NDArrayComplex]:
    """Construct basis functions for ``irrep`` by linear combinations of basis functions of ``representation``.

    Parameters
    ----------
    representation: array, (order, dim, dim)
    irrep: array, (order, dim_irrep, dim_irrep)
        Unitary (projective) irrep with factor system s.t. :math:`\\mu(E, E) = 1`.
    atol: float, default=1e-5
        Absolute tolerance to compare basis vectors

    Returns
    -------
    basis: list of array with (irrep_dim, dim)
        Each basis vectors are orthonormal.
    """
    order = irrep.shape[0]
    dim_irrep = irrep.shape[1]
    dim = representation.shape[1]
    if representation.shape != (order, dim, dim) or irrep.shape != (order, dim_irrep, dim_irrep):
        raise ValueError("Given representation and irrep do not have consistent dimensions.")

    # Pre-compute number of independent basis vectors
    character_irrep = get_character(irrep)
    character = get_character(representation)
    num_basis = np.sum(np.conj(character_irrep) * character) / order
    num_basis = np.real(np.around(num_basis)).astype(int)
    if num_basis == 0:
        return []

    count = 0
    basis: list[NDArrayComplex] = []
    for n in range(dim):
        # Initial vector to be applied projection operator
        phi = np.zeros(dim, dtype=np.complex_)
        phi[n] = 1.0

        for j in range(dim_irrep):
            # basis_nj[i, :] is the i-th basis vector forms given irrep (i = 0, ... dim_irrep-1)
            # These basis vectors are mutually orthogonal by construction!
            basis_nj = (
                dim_irrep
                / order
                * np.einsum(
                    "ki,km->im",
                    np.conj(irrep[:, :, j]),
                    representation[:, :, n],
                    optimize="greedy",
                )
            )

            if np.allclose(basis_nj, 0, atol=atol):
                continue

            # Them, normalize basis vectors s.t. they are orthonormal.
            basis_nj /= np.linalg.norm(basis_nj, axis=1)[:, None]

            # Check if linearly independent with other basis vectors
            # Two subspaces spanned by orthonormal basis vectors V1 and V2 are the same if and only if
            # triangular matrices R1 and R2 in QR decomposition of V1 and V2 are the same.
            if any([_is_same_subspace(basis_nj, other) for other in basis]):
                continue

            basis.append(basis_nj)
            count += 1
            if count == num_basis:
                break

    if count != num_basis:
        warn(
            f"Inconsistent number of independent basis vectors: expect={num_basis}, actual={count}"
        )

    return basis


def _is_same_subspace(
    basis1: NDArrayComplex,
    basis2: NDArrayComplex,
    atol: float = 1e-8,
) -> bool:
    """
    Parameters
    ----------
    basis1: array, (dim_irrep, dim)
    basis2: array, (dim_irrep, dim)
    """
    _, r1 = np.linalg.qr(basis1.T)
    _, r2 = np.linalg.qr(basis2.T)

    # If basis1 and basis2 are equivalent, r1 and r2 are the same up to U(1) phase
    for i in range(r1.shape[0]):
        if np.isclose(r2[i, i], 0, atol=atol):
            continue
        phase = r1[i, i] / r2[i, i]
        break

    return np.allclose(r1, r2 * phase, atol=atol)


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
    atol: float = 1e-8,
) -> bool:
    for i, ri in enumerate(rep):
        for j, rj in enumerate(rep):
            actual = ri @ rj
            expect = rep[table[i, j]] * factor_system[i, j]
            if not np.allclose(actual, expect, rtol=rtol, atol=atol):
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
    indicator = np.einsum("kij,kji->", irrep, irrep, optimize="greedy") / order
    indicator = int(np.around(np.real(indicator)))

    if indicator > 1:
        raise ValueError(f"Given representation is not irreducible: indicator={indicator}")

    return indicator


def check_spacegroup_representation(
    little_rotations: NDArrayInt,
    little_translations: NDArrayFloat,
    kpoint: NDArrayFloat,
    rep: NDArrayComplex,
    rtol: float = 1e-5,
) -> bool:
    """Check definition of representation. This function works for primitive and conventional cell."""
    little_rotations_int = [ndarray2d_to_integer_tuple(rotation) for rotation in little_rotations]

    # Check if ``rep`` preserves multiplication
    for r1, t1, m1 in zip(little_rotations, little_translations, rep):
        for r2, t2, m2 in zip(little_rotations, little_translations, rep):
            r12 = r1 @ r2
            t12 = r1 @ t2 + t1
            idx = little_rotations_int.index(ndarray2d_to_integer_tuple(r12))
            # little_translations[idx] may differ from t12 by lattice translation.
            m12 = rep[idx] * np.exp(-2j * np.pi * np.dot(kpoint, t12 - little_translations[idx]))

            if not np.allclose(m12, m1 @ m2, rtol=rtol):
                return False

    return True
