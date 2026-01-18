from __future__ import annotations

from itertools import product
from warnings import warn

import numpy as np

from spgrep._constants import MAX_NUM_RANDOM_GENERATIONS, RTOL
from spgrep.utils import NDArrayComplex, NDArrayFloat

from .irreps import is_equivalent_irrep


def enumerate_unitary_irreps_from_regular_representation(
    reg: NDArrayComplex,
    rtol: float = RTOL,
    max_num_random_generations: int = MAX_NUM_RANDOM_GENERATIONS,
) -> list[NDArrayComplex]:
    """Decompose given (projective) regular representation and obtain all unitary Irreps.

    Parameters
    ----------
    reg: array, (order, order, order)
        (Projective) Regular representation. reg[k] is a representation matrix for the k-th operation.
    rtol: float
        Relative tolerance to distinguish difference eigenvalues
    max_num_random_generations: int
        Maximum number of trials to generate random matrix

    Returns
    -------
    irreps: list of unitary Irreps with (order, dim, dim)
    """
    n = reg.shape[0]

    # For (m, i), reg[m, i, :] has only one nonzero entry.
    # To reduce computational time, suppress reg to only nonzero elements
    reg_nonzero = np.zeros((n, n), dtype=np.complex128)
    lookup = np.zeros((n, n), dtype=int)
    for m, i in product(range(n), repeat=2):
        idx = np.nonzero(reg[m, i, :])[0].item()
        reg_nonzero[m, i] = reg[m, i, idx]
        lookup[m, i] = idx

    rng = np.random.default_rng(seed=0)
    for _ in range(max_num_random_generations):
        # Randomly generate Hermite matrix
        hermite_random = rng.random((n, n)) + rng.random((n, n)) * 1j
        hermite_random += np.conj(hermite_random.T)

        hermite_random_reordered = np.zeros((n, n, n), dtype=np.complex128)
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
        irreps = _get_irreps_from_matrix(reg, matrix, rtol=rtol)

        if np.sum([irrep.shape[1] ** 2 for irrep in irreps]) == n:
            return irreps

    warn("Failed to search all irreps. Try increasing max_num_random_generations.")
    return []


def decompose_representation(
    representation: NDArrayComplex,
    rtol: float = RTOL,
    max_num_random_generations: int = MAX_NUM_RANDOM_GENERATIONS,
) -> list[NDArrayComplex]:
    """Decompose given (projective) representation into all unitary irreps.

    Parameters
    ----------
    representation: array, (order, dim0, dim0)
        (Projective) representation. representation[k] is a representation matrix for the k-th operation.
    rtol: float
        Relative tolerance to distinguish difference eigenvalues
    max_num_random_generations: int
        Maximum number of trials to generate random matrix

    Returns
    -------
    irreps: list of unitary Irreps with (order, dim, dim)
    """
    dim0 = representation.shape[1]

    rng = np.random.default_rng(seed=0)
    for _ in range(max_num_random_generations):
        # Randomly generate Hermite matrix
        hermite_random = rng.random((dim0, dim0)) + rng.random((dim0, dim0)) * 1j
        hermite_random += np.conj(hermite_random.T)

        # Construct matrix which commute with regular representation
        matrix = np.einsum(
            "mik,kl,mjl->ij",
            representation,
            hermite_random,
            np.conj(representation),
            optimize="greedy",
        )

        # Decompose to subspaces corresponding to Irreps
        irreps = _get_irreps_from_matrix(representation, matrix, rtol=rtol)

        if np.sum([irrep.shape[1] for irrep in irreps]) == dim0:
            return irreps

    warn("Failed to search all irreps. Try increasing max_num_random_generations.")
    return []


def _get_irreps_from_matrix(
    reg: NDArrayComplex, matrix: NDArrayComplex, rtol: float = RTOL
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
        # character = np.einsum("li,klm,mi->k", np.conj(transformation), reg, transformation)
        #           = np.einsum("klm,ml->k", reg, proj)
        proj = transformation @ np.conj(transformation.T)
        character = np.einsum("klm,ml->k", reg, proj, optimize="greedy")

        # Check if this is really irrep by character
        if not is_equivalent_irrep(character, character):
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

        # irrep = np.einsum("li,klm,mj->kij", np.conj(transformation), reg, transformation)
        tmp = reg @ transformation  # (order, order, dim)
        irrep = np.einsum("li,klj->kij", np.conj(transformation), tmp, optimize="greedy")
        irreps.append(irrep)
        characters.append(character)

    # sort Irreps by (dim, minus of sum of characters)
    argidx = sorted(range(len(irreps)), key=lambda i: (irreps[i].shape[1], -np.sum(characters[i])))
    sorted_irreps = [irreps[i] for i in argidx]
    return sorted_irreps
