from __future__ import annotations

from itertools import product
from warnings import warn

import numpy as np

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
    reg = np.zeros((n, n, n), dtype=int)
    rotations_list = [ndarray2d_to_integer_tuple(r) for r in rotations]

    for k, gk in enumerate(rotations):
        for j, gj in enumerate(rotations):
            gkj = np.dot(gk, gj)
            try:
                i = rotations_list.index(ndarray2d_to_integer_tuple(gkj))
            except ValueError:
                raise ValueError("Given matrices should form group.")

            reg[k, i, j] = 1

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
    reg = np.zeros((n, n, n), dtype=np.complex_)
    rotations_list = [ndarray2d_to_integer_tuple(r) for r in rotations]

    for k, gk in enumerate(rotations):
        for j, gj in enumerate(rotations):
            gkj = np.dot(gk, gj)
            try:
                i = rotations_list.index(ndarray2d_to_integer_tuple(gkj))
            except ValueError:
                raise ValueError("Given matrices should form group.")

            reg[k, i, j] = factor_system[k, j]

    return reg


def get_irreps(
    reg: NDArrayComplex,
    rtol: float = 1e-5,
    max_num_random_generations: int = 4,
) -> list[NDArrayComplex]:
    """Decompose given (projective) regular representation and obtain all unitary Irreps.

    Parameters
    ----------
    reg: array, (order, order, order)
        (Projective) Regular representation. reg[k] is a representation matrix for the k-th operation.
    rtol: float
        Relative tolerance to distinguish difference eigenvalues
    max_num_random_generations: int
        Maximal number of trials to generate random matrix

    Returns
    -------
    irreps: list of unitary Irreps with (order, dim, dim)
    """
    n = reg.shape[0]

    # For (m, i), reg[m, i, :] has only one nonzero entry.
    # To reduce computational time, suppress reg to only nonzero elements
    reg_nonzero = np.zeros((n, n), dtype=np.complex_)
    lookup = np.zeros((n, n), dtype=int)
    for m, i in product(range(n), repeat=2):
        idx = np.nonzero(reg[m, i, :])[0]
        reg_nonzero[m, i] = reg[m, i, idx]
        lookup[m, i] = idx

    rng = np.random.default_rng(seed=0)
    for _ in range(max_num_random_generations):
        # Randomly generate Hermite matrix
        hermite_random = rng.random((n, n)) + rng.random((n, n)) * 1j
        hermite_random += np.conj(hermite_random.T)

        hermite_random_reordered = np.zeros((n, n, n), dtype=np.complex_)
        meshi, meshj = np.meshgrid(range(n), range(n))
        # hermite_random_reordered[m, i, j] = hermite_random[lookup[m, i], lookup[m, j]]
        for m in range(n):
            hermite_random_reordered[m] = hermite_random[lookup[m, meshi], lookup[m, meshj]]

        # Construct matrix which commute with regular representation
        # Equivalent to np.einsum("mik,kl,mjl->ij", reg, hermite_random, np.conj(reg)),
        # but einsum version takes O(n^5), whereas this implementation takes O(n^3).
        # Broadcast to "mij" and sum over "m"
        matrix = np.sum(
            reg_nonzero[:, :, None] * hermite_random_reordered * np.conj(reg_nonzero[:, None, :]),
            axis=0,
        )

        # Decompose to subspaces corresponding to Irreps
        irreps = _get_irreps_from_matrix(reg, matrix, rtol)

        if np.sum([irrep.shape[1] ** 2 for irrep in irreps]) == n:
            return irreps

    warn("Failed to search all irreps. Try increasing max_num_random_generations.")
    return []


def _get_irreps_from_matrix(
    reg: NDArrayComplex, matrix: NDArrayComplex, rtol: float
) -> list[NDArrayComplex]:
    # eigvecs[:, i] is the normalized eigenvector to eigvals[i]
    eigvals, eigvecs = np.linalg.eigh(matrix)

    # Group by eigenvectors by eigenvalues
    eigenspaces: list[(float, list[NDArrayFloat])] = []  # type: ignore
    for eigval, eigvec in zip(eigvals, eigvecs.T):
        is_new_space = True
        for i, (eigval2, _) in enumerate(eigenspaces):
            if np.isclose(eigval, eigval2, rtol=rtol):
                eigenspaces[i][1].append(eigvec)
                is_new_space = False
                break
        if is_new_space:
            eigenspaces.append((eigval, [eigvec]))

    order = reg.shape[0]
    irreps: list[NDArrayComplex] = []
    characters: list[NDArrayFloat] = []
    for eigval, list_eigvecs in eigenspaces:
        # QR decomposition of column-wise vectors gives Gram-Schmidt orthonormalized vectors in column wise.
        transformation = np.linalg.qr(np.transpose(list_eigvecs))[0]

        # Compute character before irrep to avoid calculating duplicated irreps
        character = np.einsum(
            "li,klm,mi->k", np.conj(transformation), reg, transformation, optimize="greedy"
        )
        # Check if this is really irrep by character
        if not np.isclose(np.sum(np.conj(character) * character), order):
            continue

        # Multi-dimensional irreps appeared several times in `eigenspaces`.
        # Therefore, we pick one of them by checking character of each irrep.
        is_unique = True
        for character2 in characters:
            product = np.around(np.sum(np.conj(character) * character2))
            if np.isclose(product, order, rtol=rtol):
                is_unique = False
                break
        if not is_unique:
            continue

        irrep = np.einsum(
            "li,klm,mj->kij", np.conj(transformation), reg, transformation, optimize="greedy"
        )
        irreps.append(irrep)
        characters.append(character)

    # sort Irreps by (dim, minus of sum of characters)
    argidx = sorted(range(len(irreps)), key=lambda i: (irreps[i].shape[1], -np.sum(characters[i])))
    sorted_irreps = [irreps[i] for i in argidx]
    return sorted_irreps


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


def frobenius_schur_indicator(irrep: NDArrayComplex) -> int:
    """Inspect given irrep is real, pseudo-real, or not unitary equivalent.

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
    indicator = np.sum(irrep * irrep.T) / order
    indicator = int(np.around(np.real(indicator)))

    if indicator > 1:
        raise ValueError(f"Given representation is not irreducible: indicator={indicator}")

    return indicator
