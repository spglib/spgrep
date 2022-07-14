from __future__ import annotations

from itertools import product
from warnings import warn

import numpy as np

from spgrep.group import get_identity_index, get_inverse_index, get_order
from spgrep.representation import get_character, get_intertwiner
from spgrep.utils import NDArrayComplex, NDArrayFloat, NDArrayInt, is_prime, nroot


def get_irreps_from_regular(
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


def decompose_representation(
    representation: NDArrayComplex,
    rtol: float = 1e-5,
    max_num_random_generations: int = 4,
) -> list[NDArrayComplex]:
    """Decompose given (projective) representation into all unitary irreps.

    Parameters
    ----------
    representation: array, (order, dim0, dim0)
        (Projective) representation. representation[k] is a representation matrix for the k-th operation.
    rtol: float
        Relative tolerance to distinguish difference eigenvalues
    max_num_random_generations: int
        Maximal number of trials to generate random matrix

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
            "mik,kl,mjl->ij", representation, hermite_random, np.conj(representation)
        )

        # Decompose to subspaces corresponding to Irreps
        irreps = _get_irreps_from_matrix(representation, matrix, rtol)

        if np.sum([irrep.shape[1] for irrep in irreps]) == dim0:
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

        irrep = np.einsum(
            "li,klm,mj->kij", np.conj(transformation), reg, transformation, optimize="greedy"
        )
        irreps.append(irrep)
        characters.append(character)

    # sort Irreps by (dim, minus of sum of characters)
    argidx = sorted(range(len(irreps)), key=lambda i: (irreps[i].shape[1], -np.sum(characters[i])))
    sorted_irreps = [irreps[i] for i in argidx]
    return sorted_irreps


def get_irreps_from_solvable_group_chain(
    table: NDArrayInt,
    factor_system: NDArrayComplex,
    solvable_chain_generators: list[int],
    rtol: float = 1e-5,
    max_num_random_generations: int = 4,
):
    """Calculate symmetrized irreps from given chain of solvable group.

    Parameters
    ----------
    table: array, (order, order)
        Cayley table
    factor_system: array, (order, order)
    solvable_group_chain: list of single generator of coset

    rtol: float
        Relative tolerance to distinguish difference eigenvalues
    max_num_random_generations: int
        Maximal number of trials to generate random matrix

    Returns
    -------
    irreps: list of unitary projective irrep with (order, dim, dim)
    """
    identity = get_identity_index(table)
    group = [identity]  # int -> GroupIdx
    irreps = [np.ones((1, 1, 1), dtype=np.complex_)]
    for r in solvable_chain_generators:
        p = get_order(table, r)  # Should be prime number
        if not is_prime(p):
            warn("Order of generators should be prime number.")
            return []

        # Power of `r`, rm[m] = r^m
        # Power of inverse of `coset_generator`, rminv[m] = r^-m
        rm = [identity]
        rinvm = [identity]
        rinv = get_inverse_index(table, r)
        for m in range(1, p):
            rm.append(table[rm[m - 1], r])
            rinvm.append(table[rinvm[m - 1], rinv])

        # Extend group by generator `r`
        subgroup = group[:]
        group = []
        for m in range(p):
            group.extend([table[rm[m], s] for s in subgroup])
        group.sort()

        subgroup_remapping = {}  # GroupIdx -> int for `subgroup`
        for i, si in enumerate(subgroup):
            subgroup_remapping[si] = i
        group_remapping = {}  # GroupIdx -> int for `group`
        for i, gi in enumerate(group):
            group_remapping[gi] = i

        # Consider induced representation and their decomposition
        next_sub_irreps = []
        for sub_irrep in irreps:
            dim = sub_irrep.shape[1]

            # Conjugated irreps with `sub_irrep`
            conj_sub_irreps = []
            for j in range(p):
                conj_sub_irrep = []
                for s in subgroup:
                    sj = table[rinvm[j], table[s, rm[j]]]
                    conj_sub_irrep.append(
                        factor_system[s, rm[j]]
                        / factor_system[rm[j], sj]
                        * sub_irrep[subgroup_remapping[sj]]
                    )
                conj_sub_irreps.append(np.array(conj_sub_irrep))

            # Check conjugated irreps are mutually equivalent or not, and construct induced representation
            conj_characters = [get_character(conj_sub_irrep) for conj_sub_irrep in conj_sub_irreps]
            if is_equivalent_irrep(conj_characters[0], conj_characters[1]):
                # Self-conjugated case

                # Scale intertwiner s.t. intertwiner^p == identity
                intertwiner = get_intertwiner(
                    conj_sub_irreps[0], conj_sub_irreps[1], rtol, max_num_random_generations
                )
                scale = intertwiner.copy()
                for _ in range(p - 1):
                    scale = np.dot(intertwiner, scale)
                intertwiner /= scale[0, 0] ** (1 / p)

                omega = nroot(np.prod([factor_system[r, rm[m]] for m in range(1, p)]), p)
                for q in range(p):
                    omegaq = omega * np.exp(2j * np.pi * q / p)
                    delta_r = intertwiner / omegaq  # Rep. matrix for r
                    delta_rm = [
                        np.eye(intertwiner.shape[0], dtype=np.complex_)
                    ]  # delta_rm[m] is rep. matrix for r^m
                    for m in range(1, p):
                        delta_rm.append(factor_system[r, rm[m]] * delta_r @ delta_rm[m - 1])

                    next_irrep = np.zeros((len(group), dim, dim), dtype=np.complex_)
                    for m in range(p):
                        for s in subgroup:
                            idx = table[rm[m], s]
                            next_irrep[group_remapping[idx]] = (
                                factor_system[rm[m], s]
                                * delta_rm[m]
                                @ sub_irrep[subgroup_remapping[s]]
                            )
                    next_sub_irreps.append(next_irrep)
            else:
                # Mutually inequivalent
                next_irrep = np.zeros((len(group), dim * p, dim * p), dtype=np.complex_)
                for m in range(p):
                    for s in subgroup:
                        idx = table[rm[m], s]
                        for j in range(p):
                            i = (j + m) % p
                            sj = table[rinvm[j], table[s, rm[j]]]
                            next_irrep[
                                group_remapping[idx],
                                i * dim : (i + 1) * dim,
                                j * dim : (j + 1) * dim,
                            ] = (
                                factor_system[idx, rm[j]]
                                / factor_system[rm[i], sj]
                                * sub_irrep[subgroup_remapping[sj]]
                            )
                next_sub_irreps.append(next_irrep)

        # Unique irreps so far
        irreps.clear()
        sub_characters = []  # type: ignore
        for sub_irrep in next_sub_irreps:
            # Skip duplicated irrep
            character = get_character(sub_irrep)
            if any([is_equivalent_irrep(character, c) for c in sub_characters]):
                continue

            irreps.append(sub_irrep)
            sub_characters.append(character)

    if group != list(range(table.shape[0])):
        warn("Generators are not sufficient to traverse group.")
        return []

    return irreps


def is_equivalent_irrep(character1: NDArrayComplex, character2: NDArrayComplex) -> bool:
    """Return true if two irreps are equivalent."""
    order = character1.shape[0]
    if np.around(np.sum(np.conj(character1) * character2)) == order:
        return True
    else:
        return False
