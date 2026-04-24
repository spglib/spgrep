from __future__ import annotations

from typing import Literal
from warnings import warn

import numpy as np

from spgrep._constants import ATOL, MAX_NUM_RANDOM_GENERATIONS, RTOL
from spgrep.rep.enumerate import enumerate_unitary_irreps_from_regular_representation
from spgrep.rep.group import get_identity_index, get_inverse_index, get_order
from spgrep.rep.irreps import frobenius_schur_indicator, is_equivalent_irrep
from spgrep.rep.representation import get_character, get_intertwiner
from spgrep.utils import NDArrayBool, NDArrayComplex, NDArrayFloat, NDArrayInt, nroot

from .group import (
    get_cayley_table,
    get_factor_system_from_little_group,
)
from .pir import _make_physically_irreps_from_complex, purify_real_irrep_value
from .pointgroup import get_pointgroup_chain_generators
from .representation import get_projective_regular_representation


def enumerate_small_representations(
    little_rotations: NDArrayInt,
    little_translations: NDArrayFloat,
    kpoint: NDArrayFloat,
    real: bool = False,
    method: Literal["Neto", "random"] = "Neto",
    rtol: float = RTOL,
    atol: float = ATOL,
    max_num_random_generations: int = MAX_NUM_RANDOM_GENERATIONS,
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

    # Identify operations that map k -> -k (k-flipping) within the little group.
    # See reality.md: when 2k not equiv 0, the physically irreducible
    # representation is defined on the little group of +/- k. Callers that want
    # that output must pass G^{k,k-bar} (see ``get_little_group_of_pm_k``); then
    # some operations will flip k.
    flip_k = np.zeros(len(little_rotations), dtype=np.bool_)
    for i, rotation in enumerate(little_rotations):
        residual = rotation.T @ kpoint - kpoint
        residual = residual - np.rint(residual)
        if not np.allclose(residual, 0, atol=atol):
            flip_k[i] = True

    if np.any(flip_k):
        return _enumerate_pir_on_pm_k(
            little_rotations,
            little_translations,
            kpoint,
            flip_k,
            factor_system,
            method=method,
            rtol=rtol,
            atol=atol,
            max_num_random_generations=max_num_random_generations,
        )

    # Every element of the given little group stabilizes k (standard G^k path).
    # For 2k equiv 0 this is the full story. For 2k not equiv 0 with a
    # G^k-only input, we fall back to the historical behavior of forcing the
    # Re/Im block realification on the G^k irreps.
    return _make_physically_irreps_from_complex(
        irreps,
        kpoint=kpoint,
        atol=atol,
        max_num_random_generations=max_num_random_generations,
    )


def _enumerate_pir_on_pm_k(
    little_rotations: NDArrayInt,
    little_translations: NDArrayFloat,
    kpoint: NDArrayFloat,
    flip_k: NDArrayBool,
    factor_system: NDArrayComplex,
    method: Literal["Neto", "random"] = "Neto",
    rtol: float = RTOL,
    atol: float = ATOL,
    max_num_random_generations: int = MAX_NUM_RANDOM_GENERATIONS,
) -> tuple[list[NDArrayFloat], list[int]]:
    r"""Physically irreducible representations on the little group of +/- k.

    Implements the construction in :ref:`physically_irreps` for the case
    :math:`2\mathbf{k} \not\equiv \mathbf{0}`. Enumerates complex small irreps
    of :math:`\mathcal{G}^{\mathbf{k}}`, induces each one to
    :math:`\mathcal{G}^{\mathbf{k}\bar{\mathbf{k}}}` (obtaining a ``2d``-dim
    complex matrix representation), and transforms it into the
    :math:`\mathbf{k}` / :math:`-\mathbf{k}`-symmetrized real basis where all
    matrix entries become real. Finally decomposes the real representation
    into real-irreducible blocks and returns those.
    """
    xsg_indices = np.where(~flip_k)[0]
    a_indices = np.where(flip_k)[0]
    assert len(a_indices) > 0
    a0_idx = int(a_indices[0])

    table = get_cayley_table(little_rotations)
    inv_a0_idx = get_inverse_index(table, a0_idx)

    xsg_indices_list = [int(i) for i in xsg_indices]
    xsg_indices_mapping: dict[int, int] = {idx: i for i, idx in enumerate(xsg_indices_list)}
    xsg_rotations = little_rotations[xsg_indices_list]
    xsg_translations = little_translations[xsg_indices_list]
    xsg_factor_system = factor_system[np.ix_(xsg_indices_list, xsg_indices_list)]

    # Complex small irreps on G^k.
    xsg_cogroup_irreps, _ = enumerate_unitary_irreps(
        rotations=xsg_rotations,
        factor_system=xsg_factor_system,
        real=False,
        method=method,
        rtol=rtol,
        atol=atol,
        max_num_random_generations=max_num_random_generations,
    )
    xsg_phases = np.array([np.exp(-2j * np.pi * np.dot(kpoint, t)) for t in xsg_translations])
    xsg_small_irreps = [rep * xsg_phases[:, None, None] for rep in xsg_cogroup_irreps]

    real_irreps: list[NDArrayFloat] = []
    indicators: list[int] = []
    seen_characters: list[NDArrayComplex] = []

    for small_irrep in xsg_small_irreps:
        # Build the induced representation Ind_{G^k}^{G^{k,k-bar}}(Gamma).
        induced = _build_induced_representation(
            little_rotations=little_rotations,
            little_translations=little_translations,
            flip_k=flip_k,
            small_irrep=small_irrep,
            table=table,
            xsg_indices=xsg_indices_list,
            xsg_indices_mapping=xsg_indices_mapping,
            a0_idx=a0_idx,
            inv_a0_idx=inv_a0_idx,
        )

        # De-duplicate: two G^k irreps that are a0-conjugates of each other
        # give rise to the same induced rep on G^{k,k-bar}.
        character = get_character(induced)
        if any(is_equivalent_irrep(character, c) for c in seen_characters):
            continue
        seen_characters.append(character)

        # At generic k (where Ind is irreducible on G^{k,k-bar}) the induced
        # rep is a genuine complex irrep. Compute the Frobenius-Schur indicator
        # on the full space-group G^{k,k-bar} (accounting for the translation
        # subgroup) and realify using the translation-averaged intertwiner
        # (:ref:`physically_irreps`). When the indicator is +1 we can find a
        # same-dimensional real form; otherwise we Re/Im-double the complex
        # representation.
        indicator = _fs_indicator_on_pm_k(induced, little_rotations, kpoint, atol=atol)
        try:
            real_irrep = _realify_small_rep_on_pm_k(
                induced,
                little_rotations,
                kpoint,
                indicator,
                atol=atol,
                max_num_random_generations=max_num_random_generations,
            )
        except (AssertionError, RuntimeError):
            # Fallback: Re/Im block doubling always produces a valid real
            # matrix rep of G^{k,k-bar}, at the cost of doubling the dimension.
            real_irrep = _realify_via_reim_block(induced)
            indicator = 0
        real_irrep = purify_real_irrep_value(real_irrep, atol=atol)
        real_irreps.append(real_irrep)
        indicators.append(indicator)

    return real_irreps, indicators


def _fs_indicator_on_pm_k(
    induced: NDArrayComplex,
    little_rotations: NDArrayInt,
    kpoint: NDArrayFloat,
    atol: float,
) -> int:
    r"""Frobenius-Schur indicator of a small rep on G^{k,k-bar}.

    Computes ``(1 / |cogroup|) sum_i 1[(I + S_i^T) k equiv 0] * Tr(rep(i^2))``,
    where ``i^2`` is the coset representative of ``(S_i, w_i)^2``. The ``1[...]``
    factor comes from the translation sum (see :ref:`physically_irreps`).

    ``induced`` is indexed on the cogroup coset representatives of
    ``little_rotations``; ``(S_i, w_i)^2`` always lies in one of those cosets
    since G^{k,k-bar} is closed.
    """
    table = get_cayley_table(little_rotations)
    order = len(little_rotations)
    total = 0.0 + 0j
    for i in range(order):
        residual = (np.eye(3) + little_rotations[i].T) @ kpoint
        residual = residual - np.rint(residual)
        if not np.allclose(residual, 0, atol=atol):
            continue
        i2 = table[i, i]
        total += np.trace(induced[i2])
    return int(np.around(np.real(total / order)))


def _realify_via_reim_block(rep: NDArrayComplex) -> NDArrayFloat:
    """Re/Im block-double realification of a complex matrix rep."""
    order, dim, _ = rep.shape
    real = np.empty((order, 2 * dim, 2 * dim), dtype=np.float64)
    re = np.real(rep)
    im = np.imag(rep)
    real[:, :dim, :dim] = re
    real[:, :dim, dim:] = -im
    real[:, dim:, :dim] = im
    real[:, dim:, dim:] = re
    return real


def _realify_small_rep_on_pm_k(
    small_rep: NDArrayComplex,
    little_rotations: NDArrayInt,
    kpoint: NDArrayFloat,
    indicator: int,
    atol: float,
    max_num_random_generations: int,
) -> NDArrayFloat:
    r"""Realify a complex small representation of G^{k,k-bar} to a real matrix rep.

    Uses the translation-averaged intertwiner
    ``U = sum_{i: (I + S_i^T) k equiv 0} small_rep(i) B small_rep(i)^{*, dagger}``
    per reality.md, then applies ``T = U^{1/2}`` (``indicator == 1``) or the
    Re/Im block doubling (``indicator == 0``) to convert to a real matrix rep.
    """
    if indicator == 1:
        # Translation-averaged intertwiner per reality.md. We start from a
        # SYMMETRIC random B (so that U = sum M B M^{*†} comes out symmetric
        # for suitable parity reasons). The translation-projection sum
        # restricts to operations ``i`` with ``(I + S_i^T) k equiv 0``.
        order, dim, _ = small_rep.shape
        rng = np.random.default_rng()
        for _ in range(max_num_random_generations):
            B = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
            B = B + B.T  # enforce symmetry
            U = np.zeros((dim, dim), dtype=np.complex128)
            for i in range(order):
                residual = (np.eye(3) + little_rotations[i].T) @ kpoint
                residual = residual - np.rint(residual)
                if not np.allclose(residual, 0, atol=atol):
                    continue
                Mi = small_rep[i]
                U += Mi @ B @ np.conj(Mi).T
            if np.linalg.matrix_rank(U, tol=max(atol, 1e-6)) == dim:
                break
        else:
            raise RuntimeError("Failed to find non-degenerate intertwiner for PIR on G^{k,k-bar}.")

        if not np.allclose(U.T, U, atol=max(atol, 1e-6)):
            raise RuntimeError("Translation-averaged intertwiner is not symmetric.")

        det = np.linalg.det(U)
        U = U / det ** (1.0 / dim)

        # Diagonalize U = S^{-1} diag(omega) S and take T = S^{-1} diag(sqrt(omega)) S.
        eigvals, eigvecs = np.linalg.eig(U)
        real_eigvecs = []
        for eigvec in np.transpose(eigvecs):
            if not np.allclose(np.real(eigvec), 0, atol=atol):
                rv = np.real(eigvec)
            else:
                rv = np.imag(eigvec)
            real_eigvecs.append(rv / np.linalg.norm(rv))
        S = np.array(real_eigvecs)
        T = S.T @ np.diag([nroot(e, 2) for e in eigvals]) @ S
        assert np.allclose(T @ T, U, atol=max(atol, 1e-6))

        real_irrep = np.real(
            np.einsum("il,klm,mj->kij", np.conj(T), small_rep, T, optimize="greedy")
        )
        return real_irrep

    if indicator in (-1, 0):
        dim = small_rep.shape[1]
        order = small_rep.shape[0]
        real_irrep = np.empty((order, 2 * dim, 2 * dim), dtype=np.float64)
        re = np.real(small_rep)
        im = np.imag(small_rep)
        real_irrep[:, :dim, :dim] = re
        real_irrep[:, :dim, dim:] = im
        real_irrep[:, dim:, :dim] = -im
        real_irrep[:, dim:, dim:] = re
        return real_irrep

    raise ValueError(f"Unexpected Frobenius-Schur indicator: {indicator}")


def _build_induced_representation(
    little_rotations: NDArrayInt,
    little_translations: NDArrayFloat,
    flip_k: NDArrayBool,
    small_irrep: NDArrayComplex,
    table: NDArrayInt,
    xsg_indices: list[int],
    xsg_indices_mapping: dict[int, int],
    a0_idx: int,
    inv_a0_idx: int,
) -> NDArrayComplex:
    """Build ``Ind_{G^k}^{G^{k,k-bar}}(Gamma)`` over the little group of +/- k.

    Uses the standard coset decomposition ``G^{k,k-bar} = G^k sqcup a0 . G^k``.
    The induced rep has dimension ``2 * d`` where ``d = dim(Gamma)`` and is a
    genuine (complex) matrix representation of ``G^{k,k-bar}``.
    """
    order = little_rotations.shape[0]
    dim = small_irrep.shape[1]

    induced = np.zeros((order, 2 * dim, 2 * dim), dtype=np.complex128)
    # For each g in G^{k,k-bar}, compute Ind(g) using:
    #   Ind(g)[i, j] = Gamma(s_i^{-1} g s_j) if s_i^{-1} g s_j in G^k else 0,
    # with coset reps (s_1, s_2) = (e, a0).
    e_idx = get_identity_index(table)
    coset_reps = [e_idx, a0_idx]
    coset_inv = [e_idx, inv_a0_idx]
    for g in range(order):
        for i in range(2):
            for j in range(2):
                # Compute s_i^{-1} . g . s_j in the full little group.
                m = table[coset_inv[i], table[g, coset_reps[j]]]
                if not flip_k[m]:  # belongs to G^k
                    induced[
                        g,
                        i * dim : (i + 1) * dim,
                        j * dim : (j + 1) * dim,
                    ] = small_irrep[xsg_indices_mapping[int(m)]]
    return induced


def enumerate_unitary_irreps(
    rotations: NDArrayInt,
    factor_system: NDArrayComplex | None = None,
    real: bool = False,
    method: Literal["Neto", "random"] = "Neto",
    rtol: float = RTOL,
    atol: float = ATOL,
    max_num_random_generations: int = MAX_NUM_RANDOM_GENERATIONS,
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
        factor_system = np.ones((order, order), dtype=np.complex128)

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

    return _make_physically_irreps_from_complex(
        irreps,
        kpoint=None,
        atol=atol,
        max_num_random_generations=max_num_random_generations,
    )


def enumerate_unitary_irreps_from_solvable_group_chain(
    table: NDArrayInt,
    factor_system: NDArrayComplex,
    solvable_chain_generators: list[int],
    atol: float = ATOL,
    max_num_random_generations: int = MAX_NUM_RANDOM_GENERATIONS,
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
    irreps = [np.ones((1, 1, 1), dtype=np.complex128)]

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
                        np.eye(intertwiner.shape[0], dtype=np.complex128)
                    ]  # delta_rm[m] is rep. matrix for r^m
                    for m in range(1, p):
                        # D(r^m) = D(r) @ D(r^{m-1}) / mu(r, r^{m-1})
                        delta_rm.append(delta_r @ delta_rm[m - 1] / factor_system[r, rm[m - 1]])

                    next_irrep = np.zeros((len(group), dim, dim), dtype=np.complex128)
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
                next_irrep = np.zeros((len(group), dim * p, dim * p), dtype=np.complex128)
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


def purify_irrep_value(irrep: NDArrayComplex, atol: float = ATOL) -> NDArrayComplex:
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
