"""On-the-fly irreps generations."""
from __future__ import annotations

from itertools import product
from typing import Literal
from warnings import warn

import numpy as np

from spgrep.group import (
    get_cayley_table,
    get_factor_system_from_little_group,
    get_identity_index,
    get_inverse_index,
    get_order,
)
from spgrep.pointgroup import get_pointgroup_chain_generators
from spgrep.representation import (
    frobenius_schur_indicator,
    get_character,
    get_intertwiner,
    get_projective_regular_representation,
)
from spgrep.utils import NDArrayComplex, NDArrayFloat, NDArrayInt, nroot


def enumerate_small_representations(
    little_rotations: NDArrayInt,
    little_translations: NDArrayFloat,
    kpoint: NDArrayFloat,
    real: bool = False,
    method: Literal["Neto", "random"] = "Neto",
    rtol: float = 1e-5,
    atol: float = 1e-8,
    max_num_random_generations: int = 4,
) -> tuple[list[NDArrayComplex], list[int]] | tuple[list[NDArrayFloat], list[int]]:
    r"""Enumerate all unitary small representations of little group.

    Parameters
    ----------
    little_rotations: array, (order, 3, 3)
    little_translations: array, (order, 3)
    kpoint: array, (3, )
    real: bool, default=False
        If True, return irreps over real vector space (so called physically irreducible representations).
        For type-II and type-III cases, representation matrix for translation :math:`(\mathbf{E}, \mathbf{t})` is chosen as

        .. math::
           \begin{pmatrix}
           \cos (\mathbf{k} \cdot \mathbf{t}) \mathbf{1}_{d} & -\sin (\mathbf{k} \cdot \mathbf{t}) \mathbf{1}_{d} \\
           \sin (\mathbf{k} \cdot \mathbf{t}) \mathbf{1}_{d} & \cos (\mathbf{k} \cdot \mathbf{t}) \mathbf{1}_{d} \\
           \end{pmatrix}

        where :math:`\mathbf{k}` is `kpoint`.

    method: str, 'Neto' or 'random'
        'Neto': construct irreps from a fixed chain of subgroups of little co-group
        'random': construct irreps by numerically diagonalizing a random matrix commute with regular representation
    rtol: float
        Relative tolerance to distinguish difference eigenvalues
    atol: float
        Relative tolerance to compare
    max_num_random_generations: int
        Maximum number of trials to generate random matrix

    Returns
    -------
    irreps: list of unitary small representations (irreps of little group) with (order, dim, dim)
    indicators: list of int
        Frobenius-Schur indicator of composed irreps of each physically irreducible representation.
    """
    factor_system = get_factor_system_from_little_group(
        little_rotations, little_translations, kpoint
    )

    # Compute irreps of little co-group
    little_cogroup_irreps, _ = enumerate_unitary_irreps(
        little_rotations,
        factor_system,
        method=method,
        rtol=rtol,
        atol=atol,
        max_num_random_generations=max_num_random_generations,
    )

    # Small representations of little group
    irreps = []
    for rep in little_cogroup_irreps:
        phases = np.array(
            [
                np.exp(-2j * np.pi * np.dot(kpoint, translation))
                for translation in little_translations
            ]
        )
        irreps.append(rep * phases[:, None, None])

    if not real:
        indicators = [frobenius_schur_indicator(irrep) for irrep in irreps]
        return irreps, indicators

    # Physically irreducible representation
    conjugated_pairs = []
    visited = [False for _ in range(len(irreps))]
    characters = [get_character(irrep) for irrep in irreps]
    for i, ci in enumerate(characters):
        if visited[i]:
            continue
        visited[i] = True
        inequivalent = False
        for j, cj in enumerate(characters):
            if visited[j]:
                continue
            if is_equivalent_irrep(np.conj(ci), cj):
                conjugated_pairs.append((i, j))
                visited[j] = True
                inequivalent = True
                break
        if not inequivalent:
            conjugated_pairs.append((i, i))

    real_irreps = []
    indicators = []
    for conj_pair in conjugated_pairs:
        irrep = irreps[conj_pair[0]]

        indicator = frobenius_schur_indicator(irrep)
        # summation over translations becomes zero unless `2*kpoint equiv 0`
        two_kpoint = 2 * kpoint
        two_kpoint -= np.rint(two_kpoint)
        if not np.allclose(two_kpoint, 0):
            indicator = 0

        real_irrep = get_physically_irrep(
            irrep, indicator, atol=atol, max_num_random_generations=max_num_random_generations
        )
        real_irrep = purify_real_irrep_value(real_irrep, atol=atol)
        real_irreps.append(real_irrep)
        indicators.append(indicator)

    return real_irreps, indicators


def enumerate_unitary_irreps(
    rotations: NDArrayInt,
    factor_system: NDArrayComplex | None = None,
    real: bool = False,
    method: Literal["Neto", "random"] = "Neto",
    rtol: float = 1e-5,
    atol: float = 1e-8,
    max_num_random_generations: int = 4,
) -> tuple[list[NDArrayComplex] | list[NDArrayFloat], list[int]]:
    """Enumerate all unitary irreps with of matrix group ``rotations`` with ``factor_system``.

    Parameters
    ----------
    rotations: array, (order, 3, 3)
    factor_system: array, (order, order)
        If not specified, treat as ordinary representation.
    real: bool, default=False
        If True, return irreps over real vector space (so called physically irreducible representations)
    method: str, 'Neto' or 'random'
        'Neto': construct irreps from a fixed chain of subgroups of little co-group
        'random': construct irreps by numerically diagonalizing a random matrix commute with regular representation
    rtol: float
        Relative tolerance to distinguish difference eigenvalues
    atol: float
        Relative tolerance to compare
    max_num_random_generations: int
        Maximum number of trials to generate random matrix

    Returns
    -------
    irreps: list of unitary irreps with (order, dim, dim)
    indicators: list of int
        Frobenius-Schur indicator of composed irreps of each physically irreducible representation.
    """
    order = rotations.shape[0]
    if factor_system is None:
        factor_system = np.ones((order, order), dtype=np.complex_)

    if method == "Neto":
        table = get_cayley_table(rotations)
        solvable_chain_generators = get_pointgroup_chain_generators(rotations)
        irreps = enumerate_unitary_irreps_from_solvable_group_chain(
            table,
            factor_system,
            solvable_chain_generators,
            atol=atol,
            max_num_random_generations=max_num_random_generations,
        )
    elif method == "random":
        reg = get_projective_regular_representation(rotations, factor_system)
        irreps = enumerate_unitary_irreps_from_regular_representation(
            reg, rtol=rtol, max_num_random_generations=max_num_random_generations
        )
    else:
        raise ValueError(f"Unknown method to compute irreps: {method}")

    # Purify values of `irreps`.
    for irrep in irreps:
        irrep = purify_irrep_value(irrep, atol=atol)

    if not real:
        indicators = [frobenius_schur_indicator(irrep) for irrep in irreps]
        return irreps, indicators

    # Physically irreducible representation
    conjugated_pairs = []
    visited = [False for _ in range(len(irreps))]
    characters = [get_character(irrep) for irrep in irreps]
    for i, ci in enumerate(characters):
        if visited[i]:
            continue
        visited[i] = True
        inequivalent = False
        for j, cj in enumerate(characters):
            if visited[j]:
                continue
            if is_equivalent_irrep(np.conj(ci), cj):
                conjugated_pairs.append((i, j))
                visited[j] = True
                inequivalent = True
                break
        if not inequivalent:
            conjugated_pairs.append((i, i))

    real_irreps = []
    indicators = []
    for conj_pair in conjugated_pairs:
        irrep = irreps[conj_pair[0]]
        indicator = frobenius_schur_indicator(irrep)
        real_irrep = get_physically_irrep(
            irrep, indicator, atol=atol, max_num_random_generations=max_num_random_generations
        )
        real_irrep = purify_real_irrep_value(real_irrep, atol=atol)
        real_irreps.append(real_irrep)
        indicators.append(indicator)

    return real_irreps, indicators


def enumerate_unitary_irreps_from_regular_representation(
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
        Maximum number of trials to generate random matrix

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
        irreps = _get_irreps_from_matrix(reg, matrix, rtol=rtol)

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
    reg: NDArrayComplex, matrix: NDArrayComplex, rtol: float = 1e-5
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


def enumerate_unitary_irreps_from_solvable_group_chain(
    table: NDArrayInt,
    factor_system: NDArrayComplex,
    solvable_chain_generators: list[int],
    atol: float = 1e-8,
    max_num_random_generations: int = 4,
):
    r"""Calculate symmetrized irreps from given chain of solvable group.

    Parameters
    ----------
    table: array, (order, order)
        Cayley table
    factor_system: array, (order, order)
    solvable_group_chain: list of single generator of coset
        Let :math:`G_{0} := G` and :math:`G_{i} := G_{i-1} / \langle` ``solvable_chain_generators[i]`` :math:`\rangle` (i = 0, 1, ...).
        Then, :math:`G_{i}` is normal subgroup of :math:`G_{i-1}` and factor group :math:`G_{i-1}/G_{i}` is Abelian.

    atol: float
        Absolute tolerance to distinguish difference eigenvalues
    max_num_random_generations: int
        Maximum number of trials to generate random matrix

    Returns
    -------
    irreps: list of unitary projective irrep with (order, dim, dim)
    """
    identity = get_identity_index(table)
    group = [identity]  # int -> GroupIdx
    irreps = [np.ones((1, 1, 1), dtype=np.complex_)]

    # Extend subgroups from identity to whole
    for r in solvable_chain_generators[::-1]:
        # Should be power of prime number
        p = get_order(table, r)

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
        group = sorted(list(set(group)))

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
                    conj_sub_irreps[0],
                    conj_sub_irreps[1],
                    atol=atol,
                    max_num_random_generations=max_num_random_generations,
                )
                scale = intertwiner.copy()
                for _ in range(p - 1):
                    scale = np.dot(intertwiner, scale)
                intertwiner /= scale[0, 0] ** (1 / p)

                omega = 1 / nroot(np.prod([factor_system[r, rm[m]] for m in range(1, p)]), p)
                for q in range(p):
                    omegaq = omega * np.exp(2j * np.pi * q / p)
                    delta_r = intertwiner / omegaq  # Rep. matrix for r
                    delta_rm = [
                        np.eye(intertwiner.shape[0], dtype=np.complex_)
                    ]  # delta_rm[m] is rep. matrix for r^m
                    for m in range(1, p):
                        # D(r^m) = D(r) @ D(r^{m-1}) / mu(r, r^{m-1})
                        delta_rm.append(delta_r @ delta_rm[m - 1] / factor_system[r, rm[m - 1]])

                    next_irrep = np.zeros((len(group), dim, dim), dtype=np.complex_)
                    for m in range(p):
                        for s in subgroup:
                            idx = table[rm[m], s]
                            next_irrep[group_remapping[idx]] = (
                                delta_rm[m]
                                @ sub_irrep[subgroup_remapping[s]]
                                / factor_system[rm[m], s]
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


def get_physically_irrep(
    irrep: NDArrayComplex,
    indicator: int,
    atol: float = 1e-5,
    max_num_random_generations: int = 4,
) -> NDArrayFloat:
    """Compute physically irreducible representation (over real number) from given unitary irrep over complex number.

    Parameters
    ----------
    irrep: array, (order, dim, dim)
        Unitary (projective) irrep
    indicator: int
        Frobenius-Schur indicator: -1, 0, or 1
    atol: float
        Relative tolerance to compare
    max_num_random_generations: int
        Maximum number of trials to generate random matrix

    Returns
    -------
    real_irrep: array, (order, dim2, dim2)
        Physically irreducible representation composed of `irrep` and its conjugated irrep
        When `indicator==1`, `dim2 == dim`.
        When `indicator==-1` or `indicator==0`, `dim2 == 2 * dim`.
    """
    order = irrep.shape[0]
    dim = irrep.shape[1]

    # Assume kpoint is commensurated
    if indicator == 1:
        # Intertwiner with determinant=1
        conj_irrep = np.conj(irrep)
        U = get_intertwiner(
            irrep, conj_irrep, atol=atol, max_num_random_generations=max_num_random_generations
        )

        # Take real or imaginary part of eigenvectors for new basis vectors
        eigvals, eigvecs = np.linalg.eig(U)  # eigvecs[:, i] is the i-th eigenvector
        real_eigvecs = []
        for eigvec in np.transpose(eigvecs):
            if not np.allclose(np.real(eigvec), 0, atol=atol):
                real_eigvec = np.real(eigvec)
            else:
                real_eigvec = np.imag(eigvec)
            real_eigvecs.append(real_eigvec / np.linalg.norm(real_eigvec))
        S = np.transpose(real_eigvecs)

        # Square root of intertwiner
        T = S.T @ np.diag([nroot(eigval, 2) for eigval in eigvals]) @ S

        real_irrep = np.real(np.einsum("il,klm,mj->kij", T, irrep, np.conj(T), optimize="greedy"))

    elif indicator in [-1, 0]:
        real_irrep = np.empty((order, 2 * dim, 2 * dim), dtype=np.float_)
        # [ [Re D(g),  Im D(g)]
        #   [-Im D(g), Re D(g)] ]
        real_irrep[:, :dim, :dim] = np.real(irrep)
        real_irrep[:, :dim, dim:] = np.imag(irrep)
        real_irrep[:, dim:, :dim] = -np.imag(irrep)
        real_irrep[:, dim:, dim:] = np.real(irrep)

    return real_irrep


def is_equivalent_irrep(character1: NDArrayComplex, character2: NDArrayComplex) -> bool:
    """Return true if two irreps are equivalent."""
    order = character1.shape[0]
    if np.around(np.sum(np.conj(character1) * character2)) == order:
        return True
    else:
        return False


def purify_irrep_value(irrep: NDArrayComplex, atol: float = 1e-8) -> NDArrayComplex:
    """Purify values of irreps."""
    # Each value should be 0 or exp(2 pi q / p) (p=1,2,3,4,6, q = 0,...,p-1)
    possible_values = [
        0,
        1,  # 0/1
        np.exp(1j * np.pi / 3),  # 1/3
        1j,  # 1/4
        np.exp(1j * np.pi * 2 / 3),  # 2/3
        -1,  # 1/2
        np.exp(1j * np.pi * 4 / 3),  # 4/3
        -1j,  # 3/4
        np.exp(1j * np.pi * 5 / 3),  # 5/3
    ]
    for v in possible_values:
        irrep[np.abs(irrep - v) < atol] = v
    return irrep


def purify_real_irrep_value(real_irrep: NDArrayFloat, atol: float = 1e-8) -> NDArrayFloat:
    """Purify values of physically irreducible representations."""
    values = [
        0,
        1,  # 0/1
        1 / 2,  # 1/3
        np.sqrt(3) / 2,  # 1/3
        -1 / 2,  # 2/3
        -np.sqrt(3) / 2,  # 2/3
        -1,  # 1/2
    ]
    for v in values:
        real_irrep[np.abs(real_irrep - v) < atol] = v
    return real_irrep
