"""Representation-matrix related implementations."""
from __future__ import annotations

from itertools import product
from warnings import warn

import numpy as np

from spgrep.group import get_cayley_table
from spgrep.utils import (
    NDArrayComplex,
    NDArrayFloat,
    NDArrayInt,
    contain_space,
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
        matrix = np.einsum("kil,lm,kjm->ij", rep1, random, np.conj(rep2), optimize="greedy")
        if not np.allclose(matrix, 0, atol=atol):
            # Scale such that determinant is unity
            matrix /= np.linalg.det(matrix) ** (1 / dim)
            return matrix

    warn("Failed to search all irreps. Try increasing max_num_random_generations.")
    return np.zeros((dim, dim))


def get_character(representation: NDArrayComplex) -> NDArrayComplex:
    """Calculate character of representation.

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
    atol: float = 1e-6,  # A bit large tolerance setting to handle numerical noise in `representation`
    max_num_trials: int = 10,
) -> list[NDArrayComplex]:
    r"""Construct basis functions for ``irrep`` by linear combinations of basis functions of ``representation``.

    Parameters
    ----------
    representation: array, (order, dim, dim)
    irrep: array, (order, dim_irrep, dim_irrep)
        Unitary (projective) irrep with factor system s.t. :math:`\mu(E, E) = 1`.
    atol: float, default=1e-5
        Absolute tolerance to compare basis vectors
    max_num_trials: int, default=10
        Maximum number to retry when failed to select projected basis vectors

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
    character_sum = np.sum(np.conj(character_irrep) * character)
    if not np.isclose(character_sum, np.around(character_sum), atol=atol):
        warn("Inner product of characters should return an integer.")
    num_basis = np.around(character_sum) / order
    num_basis = np.around(np.real(num_basis)).astype(int)
    if num_basis == 0:
        return []

    def _project_to_irrep(adjusted_atol):
        count = 0
        basis: list[NDArrayComplex] = []
        for n in range(dim):
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

                if np.allclose(basis_nj, 0, atol=adjusted_atol):
                    continue

                # Them, normalize basis vectors s.t. they are orthonormal.
                basis_nj /= np.linalg.norm(basis_nj, axis=1)[:, None]

                # Check if linearly independent with other basis vectors
                # Two subspaces spanned by orthonormal basis vectors V1 and V2 are the same if and only if
                # triangular matrices R1 and R2 in QR decomposition of V1 and V2 are the same.
                if len(basis) > 0 and contain_space(
                    np.concatenate(basis, axis=0), basis_nj, atol=adjusted_atol
                ):
                    # if any([contain_space(basis_nj, other, atol=adjusted_atol) for other in basis]):
                    continue

                basis.append(basis_nj)
                count += 1

        return basis, count

    adjusted_atol = atol
    for _ in range(max_num_trials):
        basis, count = _project_to_irrep(adjusted_atol)
        if count == num_basis:
            break
        elif count > num_basis:
            # Tighten tolerance to compare basis vectors
            adjusted_atol /= 2
            print(count, num_basis, adjusted_atol)
        else:
            # Loosen tolerance to compare basis vectors
            adjusted_atol *= 2

    if count > num_basis:
        warn(
            f"Inconsistent number of independent basis vectors (expect={num_basis}, actual={count})."
            "Try decreasing atol."
        )
    elif count < num_basis:
        warn(
            f"Inconsistent number of independent basis vectors (expect={num_basis}, actual={count})."
            "Try increasing atol."
        )

    return basis


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
    rtol: float = 1e-5,
    atol: float = 1e-8,
) -> bool:
    """Return true if given matrix function is a (projective) representation with given factor system."""
    order = rep.shape[0]
    if factor_system is None:
        factor_system = np.ones((order, order), dtype=np.complex_)

    for i, ri in enumerate(rep):
        for j, rj in enumerate(rep):
            actual = ri @ rj
            expect = rep[table[i, j]] * factor_system[i, j]
            if not np.allclose(actual, expect, rtol=rtol, atol=atol):
                return False

    return True


def frobenius_schur_indicator(irrep: NDArrayComplex) -> int:
    r"""Inspect given unitary (projective) irrep is real, pseudo-real, or not unitary equivalent.

    .. math::
       \mathrm{indicator} =
       \frac{1}{|G|} \sum_{ g \in G } \chi(g^{2})

    Parameters
    ----------
    irrep: array, (order, dim, dim)

    Returns
    -------
    indicator: int
        If ``indicator==1``, it is real Reps.
        If ``indicator==-1``, it is psedu-real Reps.
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
