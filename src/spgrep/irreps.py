from __future__ import annotations

from itertools import product
from warnings import warn

import numpy as np

from spgrep.group import (
    get_cayley_table,
    get_identity_index,
    get_inverse_index,
    get_order,
)
from spgrep.utils import NDArrayComplex, NDArrayFloat, NDArrayInt


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


def symmetrize_irrep(
    irrep: NDArrayComplex,
    table: NDArrayInt,
    factor_system: NDArrayComplex,
    solvable_group_chain: list[tuple[list[int], int]],
    rtol: float = 1e-5,
    max_num_random_generations: int = 4,
):
    """Standardize given unitary (projective) irrep.

    Parameters
    ----------
    irrep: array, (order, dim, dim)
    table: array, (order, order)
        Cayley table
    factor_system: array, (order, order)
    solvable_group_chain: list of (indices of subgroup, single generator of coset, index of coset decomposition)
        TODO
    rtol: float
        Relative tolerance to distinguish difference eigenvalues
    max_num_random_generations: int
        Maximal number of trials to generate random matrix

    Returns
    -------
    symmetrized_irrep: array, (order, dim, dim)
    """
    order = irrep.shape[0]
    return _symmetrize_irrep(
        irrep=irrep,
        table=table,
        factor_system=factor_system,
        group=list(range(order)),
        solvable_group_chain=solvable_group_chain,
        rtol=rtol,
        max_num_random_generations=max_num_random_generations,
    )


def _symmetrize_irrep(
    irrep: NDArrayComplex,
    table: NDArrayInt,
    factor_system: NDArrayComplex,
    group: list[int],  # int -> GroupIdx
    solvable_group_chain: list[tuple[list[int], int]],
    rtol: float = 1e-5,
    max_num_random_generations: int = 4,
):
    dim = irrep.shape[1]
    if dim == 1:
        # One-dimensional unitary irrep is unique
        return irrep

    # Let the i-th operation be ops[i] (i: GroupIdx).
    # Let i = group[si], irrep[i] is for the i-th operation.
    subgroup, r = solvable_group_chain[0]  # subgroup: int -> GroupIdx
    p = get_order(table, r)  # Should be prime number

    # GroupIdx -> int for `group`
    group_remapping = {}
    for i, gi in enumerate(group):
        group_remapping[gi] = i
    subgroup_remapped = [group_remapping[si] for si in subgroup]  # list[int]

    # GroupIdx -> int for `subgroup`
    subgroup_remapping = {}
    for i, si in enumerate(subgroup):
        subgroup_remapping[si] = i

    decomposed = decompose_representation(
        irrep[subgroup_remapped], rtol, max_num_random_generations
    )

    # Choose one of irreps decomposed from subduced representation
    if len(decomposed) == 1:
        subduced_irrep = decomposed[0]
    else:
        # TODO: triply degenerated case
        if len(decomposed) % 3 == 0:
            raise NotImplementedError

        # Sort decomposed irreps by angle of character of coset representative
        characters = [get_character(ir) for ir in decomposed]
        next_coset_generator = subgroup_remapping[solvable_group_chain[1][1]]
        angles = np.angle([character[next_coset_generator] for character in characters])

        subduced_irrep = decomposed[np.argmin(angles)]

    # Symmetrize irrep from subduced representation recursively
    if subduced_irrep.shape[1] == 1:
        # One-dimensional unitary irrep is unique
        symmetrized_subduced_irrep = subduced_irrep
    else:
        symmetrized_subduced_irrep = _symmetrize_irrep(
            subduced_irrep,  # (sub_order/gen_order, dim0, dim0)
            table=table,
            factor_system=factor_system,
            group=subgroup,
            solvable_group_chain=solvable_group_chain[1:],
            rtol=rtol,
            max_num_random_generations=max_num_random_generations,
        )

    # Power of `coset_generator`, rm[m] = r^m
    rm = [get_identity_index(table)]
    # Power of inverse of `coset_generator`, rminv[m] = r^-m
    rinv = get_inverse_index(table, r)
    rinvm = [get_identity_index(table)]
    for m in range(1, p):
        rm.append(table[rm[m - 1], r])
        rinvm.append(table[rinvm[m - 1], rinv])

    # Conjugated irreps with `symmetrized_subduced_irrep`
    conjugated_irreps = []
    for j in range(p):
        conj_irrep = []
        for sk in subgroup:
            sj = table[rinvm[j], table[sk, rm[j]]]
            conj_irrep.append(
                factor_system[sk, rm[j]]
                / factor_system[rm[j], sj]
                * symmetrized_subduced_irrep[subgroup_remapping[sj]]
            )
        conjugated_irreps.append(np.array(conj_irrep))

    # Check conjugated irreps are mutually equivalent or not, and construct induced representation
    conjugated_characters = [get_character(conj_irrep) for conj_irrep in conjugated_irreps]
    if is_equivalent_irrep(conjugated_characters[0], conjugated_characters[1]):
        # Self conjugated
        intertwiner = get_intertwiner(conjugated_irreps[0], conjugated_irreps[1])

        # Scale intertwiner s.t. intertwiner^p == identity
        scale = intertwiner.copy()
        for _ in range(p - 1):
            scale = np.dot(intertwiner, scale)
        intertwiner /= scale[0, 0] ** (1 / p)

        omega = np.prod([factor_system[r, rm[m]] for m in range(1, p)]) ** (1 / p)
        character = get_character(irrep)
        found = False
        for q in range(p):
            omegaq = omega * np.exp(2j * np.pi * q / p)
            delta_r = intertwiner / omegaq  # Rep. matrix for r
            delta_rm = [
                np.eye(intertwiner.shape[0], dtype=np.complex_)
            ]  # delta_rm[m] is rep. matrix for r^m
            for m in range(1, p):
                delta_rm.append(factor_system[r, rm[m]] * delta_r @ delta_rm[m - 1])

            symmetrized_irrep = np.zeros((len(group), dim, dim), dtype=np.complex_)
            for m in range(p):
                for sk in subgroup:
                    idx = table[rm[m], sk]
                    symmetrized_irrep[group_remapping[idx]] = (
                        factor_system[rm[m], sk]
                        * delta_rm[m]
                        @ symmetrized_subduced_irrep[subgroup_remapping[sk]]
                    )

            if is_equivalent_irrep(character, get_character(symmetrized_irrep)):
                found = True
                break
        if not found:
            raise ValueError("Unreachable!")
    else:
        # Mutually inequivalent
        symmetrized_irrep = np.zeros((len(group), dim, dim), dtype=np.complex_)
        dim_sub = dim // p
        for m in range(p):
            for sk in subgroup:
                idx = table[rm[m], sk]
                for j in range(p):
                    i = (j + m) % p
                    sj = table[rinvm[j], table[sk, rm[j]]]
                    symmetrized_irrep[
                        group_remapping[idx],
                        i * dim_sub : (i + 1) * dim_sub,
                        j * dim_sub : (j + 1) * dim_sub,
                    ] = (
                        factor_system[idx, rm[j]]
                        / factor_system[rm[i], sj]
                        * symmetrized_subduced_irrep[subgroup_remapping[sj]]
                    )

    # Symmetrize given irrep by generating induced representation
    return symmetrized_irrep


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


def is_equivalent_irrep(character1: NDArrayComplex, character2: NDArrayComplex) -> bool:
    """Return true if two irreps are equivalent."""
    order = character1.shape[0]
    if np.around(np.sum(np.conj(character1) * character2)) == order:
        return True
    else:
        return False


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
