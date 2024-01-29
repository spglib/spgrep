"""Spin co-representation."""

from __future__ import annotations

from typing import Literal

import numpy as np

from spgrep.group import (
    decompose_by_maximal_space_subgroup,
    get_cayley_table,
    get_factor_system_from_little_group,
    get_inverse_index,
)
from spgrep.irreps import enumerate_unitary_irreps, is_equivalent_irrep
from spgrep.representation import get_character, get_intertwiner
from spgrep.spinor import (
    enumerate_spinor_small_representations,
    get_spinor_unitary_rotation,
)
from spgrep.utils import NDArrayBool, NDArrayComplex, NDArrayFloat, NDArrayInt


def enumerate_spinor_small_corepresentations(
    lattice: NDArrayFloat,
    little_rotations: NDArrayInt,
    little_translations: NDArrayFloat,
    little_time_reversals: NDArrayInt,
    kpoint: NDArrayFloat,
    method: Literal["Neto", "random"] = "Neto",
    rtol: float = 1e-5,
    atol: float = 1e-8,
    max_num_random_generations: int = 4,
) -> tuple[list[NDArrayComplex], list[int], NDArrayComplex, NDArrayComplex, NDArrayBool]:
    r"""Enumerate all unitary co-irreps of little group for spinor.

    .. math::
       \overline{\mathbf{D}}^{\mathbf{k}\alpha}(\mathbf{S}_{i}1) \overline{\mathbf{D}}^{\mathbf{k}\alpha}(\mathbf{S}_{j}\theta_{j})
            = \omega(\mathbf{S}_{i}1, \mathbf{S}_{j}\theta_{j}) e^{ -i \mathbf{g}_{i} \cdot \mathbf{w}_{j} } \overline{\mathbf{D}}^{\mathbf{k}\alpha}(\mathbf{S}_{k}\theta_{k}) \\
       \overline{\mathbf{D}}^{\mathbf{k}\alpha}(\mathbf{S}_{i}1') \overline{\mathbf{D}}^{\mathbf{k}\alpha}(\mathbf{S}_{j}\theta_{j})^{\ast}
            = \omega(\mathbf{S}_{i}1', \mathbf{S}_{j}\theta_{j}) e^{ -i \mathbf{g}_{i} \cdot \mathbf{w}_{j} } \overline{\mathbf{D}}^{\mathbf{k}\alpha}(\mathbf{S}_{k}\theta_{k}) \\

    See :ref:`corep_spinor_factor_system` for spinor-derived factor system :math:`\omega`.

    Parameters
    ----------
    lattice: array, (3, 3)
        Row-wise basis vectors. ``lattice[i, :]`` is the i-th lattice vector.
    little_rotations: array[int], (order, 3, 3)
    little_translations: array, (order, 3)
    little_time_reversals: array[int], (order)
    kpoint: array, (3, )
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
    small_coreps: list of array
        Unitary small co-representations (irreps of little group), :math:`\{ \overline{\mathbf{D}}^{\mathbf{k}\alpha} \}_{\alpha}`
    indicators: list[int]
        Frobenius-Schur indicators for each small representation
    corep_spinor_factor_system: array, (order, order)
        ``factor_system[i, j]`` stands for :math:`\omega(\mathbf{S}_{i}\theta_{i}, \mathbf{S}_{j}\theta_{j})`
    unitary_rotations: array, (order, 2, 2)
        ``unitary_rotations[i]`` stands for :math:`\mathbf{U}(\mathbf{S}_{i}) \in SU(2)`.
    anti_linear: array[bool], (order, )
        If ``anti_linear[i] == True``, the ``i``-th operator is anti-linear.
    """
    order = len(little_rotations)

    if np.all(little_time_reversals == 0):
        # Type-I MSG
        irreps, spinor_factor_system, unitary_rotations = enumerate_spinor_small_representations(
            lattice=lattice,
            little_rotations=little_rotations,
            little_translations=little_translations,
            kpoint=kpoint,
            method=method,
            rtol=rtol,
            atol=atol,
            max_num_random_generations=max_num_random_generations,
        )
        anti_linear = np.zeros((order,), dtype=np.bool_)
        return irreps, [], spinor_factor_system, unitary_rotations, anti_linear

    # Coset-decompose Type-II, -III, -IV MSG
    (xsg_indices, _, a0_idx) = decompose_by_maximal_space_subgroup(  # type: ignore
        little_rotations, little_translations, little_time_reversals
    )
    xsg_order = len(xsg_indices)
    xsg_indices_mapping = {}  # xsg_indices -> [0, xsg_order)
    for i, idx in enumerate(xsg_indices):
        xsg_indices_mapping[idx] = i

    # Factor system from spinor
    corep_spinor_factor_system, unitary_rotations, anti_linear = get_corep_spinor_factor_system(
        lattice, little_rotations, little_time_reversals
    )
    # Factor system from nonsymmorphic
    nonsymmorphic_factor_system = get_factor_system_from_little_group(
        little_rotations,
        little_translations,
        kpoint,
    )
    factor_system = corep_spinor_factor_system * nonsymmorphic_factor_system

    # Compute irreps of little co-group
    xsg_factor_system = factor_system[xsg_indices][:, xsg_indices]
    little_cogroup_irreps, _ = enumerate_unitary_irreps(
        rotations=little_rotations[xsg_indices],
        factor_system=xsg_factor_system,
        real=False,  # Nonsense to consider real-value irreps
        method=method,
        rtol=rtol,
        atol=atol,
        max_num_random_generations=max_num_random_generations,
    )

    # Construct "conjugated" irreps
    conj_irreps = []
    table = get_cayley_table(little_rotations, little_time_reversals)
    characters = [get_character(irrep) for irrep in little_cogroup_irreps]
    inv_a0_idx = get_inverse_index(table, a0_idx)
    for irrep in little_cogroup_irreps:
        conj_indices = table[inv_a0_idx, table[xsg_indices, a0_idx]]  # a0^-1 * xsg_indices * a0
        conj_indices_mapping = [xsg_indices_mapping[idx] for idx in conj_indices]
        conj_irrep = (
            factor_system[xsg_indices, a0_idx][:, None, None]
            / factor_system[a0_idx, conj_indices][:, None, None]
            * np.conj(irrep[conj_indices_mapping])
        )
        conj_irreps.append(conj_irrep)

    # Pair "conjugated" irreps
    conj_pairs = []  # list of (indicator, irrep of little co-group, conjugated irrep)
    indicators = []
    visited = [False for _ in range(len(little_cogroup_irreps))]
    a0u = table[a0_idx, xsg_indices]  # a0 * u
    factor_a0u_a0u = [factor_system[idx, idx] for idx in a0u]
    a0ua0u_mapping = [xsg_indices_mapping[table[idx, idx]] for idx in a0u]  # (a0 * u)^2
    for i, (character_i, irrep_i, conj_irrep_i) in enumerate(
        zip(characters, little_cogroup_irreps, conj_irreps)
    ):
        if visited[i]:
            continue
        visited[i] = True

        # Frobenius-Schur indicator for co-representation
        indicator = np.around(
            np.real(np.sum(factor_a0u_a0u * character_i[a0ua0u_mapping]))
        ).astype(int)
        assert indicator % xsg_order == 0
        indicator = indicator // xsg_order

        conj_pairs.append((indicator, irrep_i, conj_irrep_i))
        indicators.append(indicator)
        if indicator != 0:
            continue  # indicator = 1 or -1

        # Inequivalent case
        conj_character_i = get_character(conj_irrep_i)
        found = False
        for j, character_j in enumerate(characters):
            if visited[j]:
                continue
            if is_equivalent_irrep(conj_character_i, character_j):
                found = True
                visited[j] = True
                break
        assert found

    # Construct co-representations
    small_coreps = []
    phases = np.array(
        [np.exp(-2j * np.pi * np.dot(kpoint, translation)) for translation in little_translations]
    )
    for indicator, irrep, conj_irrep in conj_pairs:
        dim = irrep.shape[1]
        if indicator == 1:
            # Unitary matrix s.t. irrep @ U = conj_irrep @ U
            U = get_intertwiner(irrep, conj_irrep, atol, max_num_random_generations)
            corep = np.zeros((order, dim, dim), dtype=np.complex_)
            corep[xsg_indices] = irrep
            corep[a0u] = (
                np.conj(factor_system[a0_idx, xsg_indices])[:, None, None] * U[None, :, :] @ irrep
            )
        elif indicator == -1:
            U = get_intertwiner(irrep, conj_irrep, atol, max_num_random_generations)
            corep = np.zeros((order, 2 * dim, 2 * dim), dtype=np.complex_)

            # [ [irrep, 0],
            #   [0, conj_irrep]]
            corep[xsg_indices, :dim, :dim] = irrep
            corep[xsg_indices, dim:, dim:] = conj_irrep

            # [ [0, -U],
            #   [U, 0] ]
            corep_a0 = np.zeros((2 * dim, 2 * dim), dtype=np.complex_)
            corep_a0[:dim, dim:] = -U
            corep_a0[dim:, :dim] = U

            corep[a0u] = (
                np.conj(factor_system[a0_idx, xsg_indices])[:, None, None]
                * corep_a0[None, :, :]
                @ corep[xsg_indices]
            )
        elif indicator == 0:
            corep = np.zeros((order, 2 * dim, 2 * dim), dtype=np.complex_)

            # [ [irrep, 0],
            #   [0, conj_irrep]]
            corep[xsg_indices, :dim, :dim] = irrep
            corep[xsg_indices, dim:, dim:] = conj_irrep

            # [ [0, omega(a0, a0) irrep[a0 * a0]],
            #   [1, 0] ]
            corep_a0 = np.zeros((2 * dim, 2 * dim), dtype=np.complex_)
            corep_a0[:dim, dim:] = (
                factor_system[a0_idx, a0_idx] * irrep[xsg_indices_mapping[table[a0_idx, a0_idx]]]
            )
            corep_a0[dim:, :dim] = np.eye(dim, dtype=np.complex_)

            corep[a0u] = (
                np.conj(factor_system[a0_idx, xsg_indices])[:, None, None]
                * corep_a0[None, :, :]
                @ corep[xsg_indices]
            )
        else:
            raise ValueError("Unreachable!")

        # Small co-representation
        small_corep = corep * phases[:, None, None]
        small_coreps.append(small_corep)

    return small_coreps, indicators, corep_spinor_factor_system, unitary_rotations, anti_linear


def get_corep_spinor_factor_system(
    lattice: NDArrayFloat,
    rotations: NDArrayInt,
    time_reversals: NDArrayInt,
) -> tuple[NDArrayComplex, NDArrayComplex, NDArrayBool]:
    r"""Calculate spin-derived factor system of spin co-representation.

    See :ref:`corep_spinor_factor_system` for spinor-derived factor system :math:`\omega`.
    An ordinary symmetry operation :math:`\mathbf{S}_{i}1` maps to unitary operator :math:`\mathbf{U}(\mathbf{S}_{i})`.
    An antisymmetry operation :math:`\mathbf{S}_{i}1'` maps to

    Parameters
    ----------
    lattice: array, (3, 3)
        Row-wise basis vectors. ``lattice[i, :]`` is the i-th lattice vector.
    rotations: array[int], (order, 3, 3)
        Rotation parts of magnetic point group, :math:`\mathbf{S}_{i}`
    time_reversals: array[int], (order, )
        Time-reversal parts of magnetic point group, :math:`\theta_{i}`

    Returns
    -------
    corep_spinor_factor_system: array, (order, order)
        ``factor_system[i, j]`` stands for :math:`\omega(\mathbf{S}_{i}\theta_{i}, \mathbf{S}_{j}\theta_{j})`
    unitary_rotations: array, (order, 2, 2)
        ``unitary_rotations[i]`` stands for :math:`\mathbf{U}(\mathbf{S}_{i}) \in SU(2)`.
    anti_linear: array[bool], (order, )
        If ``anti_linear[i] == True``, the ``i``-th operator is anti-linear.
    """
    # -i sigma_y
    time_reversal_matrix = np.array(
        [
            [0, -1],
            [1, 0],
        ],
        dtype=np.complex_,
    )

    # Assign a unitary or anti-unitary operator for each magnetic operation
    order = len(rotations)
    unitary_rotations = np.zeros((order, 2, 2), dtype=np.complex_)
    anti_linear = np.zeros((order,), dtype=np.bool_)
    for i, (rotation, time_reversal) in enumerate(zip(rotations, time_reversals)):
        unitary_rotation = get_spinor_unitary_rotation(lattice, rotation)
        if time_reversal == 1:
            # Anti-linear operator
            unitary_rotations[i] = unitary_rotation @ time_reversal_matrix
            anti_linear[i] = True
        else:
            # Linear operator
            unitary_rotations[i] = unitary_rotation
            anti_linear[i] = False

    # Factor system for spinor co-rep
    table = get_cayley_table(rotations, time_reversals)
    corep_spinor_factor_system = np.zeros((order, order), dtype=np.complex_)
    for i, (ui, ai) in enumerate(zip(unitary_rotations, anti_linear)):
        for j, uj in enumerate(unitary_rotations):
            if ai:
                uiuj = ui @ np.conj(uj)  # ui is anti-linear
            else:
                uiuj = ui @ uj

            # si @ sj = sk in O(3)
            k = table[i, j]
            # Multiplier should be -1 or 1
            for multiplier in [-1, 1]:
                if np.allclose(uiuj, multiplier * unitary_rotations[k]):
                    corep_spinor_factor_system[i, j] = multiplier
                    break

    return corep_spinor_factor_system, unitary_rotations, anti_linear
