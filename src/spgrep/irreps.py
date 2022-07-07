from __future__ import annotations

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
    """Decompose given (projective) regular representation and obtain all Irreps.

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
    irreps: list of Irreps with (order, dim, dim)
    """
    n = reg.shape[0]

    rng = np.random.default_rng(seed=0)
    for _ in range(max_num_random_generations):
        # Randomly generate Hermite matrix
        hermite_random = rng.random((n, n)) + rng.random((n, n)) * 1j
        hermite_random += np.conj(hermite_random.T)

        # Construct matrix which commute with regular representation
        matrix = np.einsum("mik,kl,mjl->ij", reg, hermite_random, np.conj(reg))

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
        # Stack eigenvectors in column wise
        transformation = np.transpose(list_eigvecs)
        irrep = np.einsum("li,klm,mj->kij", np.conj(transformation), reg, transformation)
        character = get_character(irrep)

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
        if is_unique:
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
    character = np.einsum("ijj->i", representation)
    return character
