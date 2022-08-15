"""Enumerate unitary irreps from crystal structure or symmetry operations."""
from __future__ import annotations

from typing import Literal

import numpy as np
from spglib import get_symmetry_dataset

from spgrep.group import get_little_group
from spgrep.irreps import (
    enumerate_small_representations,
    enumerate_unitary_irreps,
    purify_irrep_value,
)
from spgrep.transform import (
    get_primitive_transformation_matrix,
    transform_symmetry_and_kpoint,
    unique_primitive_symmetry,
)
from spgrep.utils import NDArrayComplex, NDArrayFloat, NDArrayInt


def get_spacegroup_irreps(
    lattice: NDArrayFloat,
    positions: NDArrayFloat,
    numbers: NDArrayInt,
    kpoint: NDArrayFloat,
    method: Literal["Neto", "random"] = "Neto",
    reciprocal_lattice: NDArrayFloat | None = None,
    symprec: float = 1e-5,
    rtol: float = 1e-5,
    atol: float = 1e-8,
    max_num_random_generations: int = 4,
) -> tuple[list[NDArrayComplex], NDArrayInt, NDArrayFloat, NDArrayInt]:
    r"""Compute all irreducible representations of space group of given structure up to unitary transformation.

    Parameters
    ----------
    lattice: array, (3, 3)
        Row-wise basis vectors. ``lattice[i, :]`` is the i-th lattice vector.
    positions: array, (num_atoms, 3)
        Fractional coordinates of sites
    numbers: array, (num_atoms, )
        Integer list specifying atomic species
    method: str, 'Neto' or 'random'
        'Neto': construct irreps from a fixed chain of subgroups of little co-group
        'random': construct irreps by numerically diagonalizing a random matrix commute with regular representation
    kpoint: array, (3, )
        Reciprocal vector with respect to ``reciprocal_lattice``
        For pure translation :math:`\mathbf{t}`, returned irrep :math:`\Gamma^{(\alpha)}` takes

        .. math::
            \Gamma^{(\alpha)}((E, \mathbf{t})) = e^{ -i\mathbf{k}\cdot\mathbf{t} } \mathbf{1}.

    reciprocal_lattice: (Optional) array, (3, 3)
        ``reciprocal_lattice[i, :]`` is the i-th basis vector of reciprocal lattice for ``kpoint`` without `2 * pi factor`.
        If not specified, ``reciprocal_lattice`` is set to ``np.linalg.inv(lattice).T``.
    symprec: float
        Parameter for searching symmetry operation in Spglib
    rtol: float
        Relative tolerance for comparing float values
    atol: float
        Absolute tolerance to distinguish difference eigenvalues
    max_num_random_generations: int
        Maximum number of trials to generate random matrix

    Returns
    -------
    irreps: list of Irreps with (little_group_order, dim, dim)
        ``irreps[alpha][i, :, :]`` is the ``alpha``-th irreducible matrix representation of ``(little_rotations[i], little_translations[i])``.
    rotations: array, (num_sym, 3, 3)
    translations: array, (num_sym, 3)
    mapping_little_group: array, (little_group_order, )
        Let ``i = mapping_little_group[idx]``.
        (rotations[i], translations[i]) belongs to the little group of given space space group and kpoint.
    """
    # Transform given `kpoint` in dual of `lattice`
    dual_lattice = np.linalg.inv(lattice).T
    if reciprocal_lattice is None:
        reciprocal_lattice = dual_lattice
    # kpoint @ reciprocal_lattice == kpoint_conv @ dual_lattice
    kpoint_conv = kpoint @ reciprocal_lattice @ np.linalg.inv(dual_lattice)

    dataset = get_symmetry_dataset(cell=(lattice, positions, numbers), symprec=symprec)
    rotations = dataset["rotations"]
    translations = dataset["translations"]

    # Transform to primitive
    to_primitive = get_primitive_transformation_matrix(dataset["hall_number"])
    prim_rotations, prim_translations, prim_kpoint = transform_symmetry_and_kpoint(
        to_primitive, rotations, translations, kpoint_conv
    )
    # mapping_to_prim: [0..num_sym) -> [0..order)
    uniq_prim_rotations, uniq_prim_translations, mapping_to_prim = unique_primitive_symmetry(
        prim_rotations, prim_translations
    )

    # mapping_prim_little_group: [0..prim_little_group_order) -> [0..order)
    prim_irreps, mapping_prim_little_group = get_spacegroup_irreps_from_primitive_symmetry(
        rotations=uniq_prim_rotations,
        translations=uniq_prim_translations,
        kpoint=prim_kpoint,
        method=method,
        rtol=rtol,
        atol=atol,
        max_num_random_generations=max_num_random_generations,
    )
    remapping_prim_little_group = {}  # [0..order) -> [0..prim_little_group_order)
    for i, idx in enumerate(mapping_prim_little_group):
        remapping_prim_little_group[idx] = i

    mapping_little_group = []  # [0..little_group_order) -> [0..num_sym)
    mapping_conv_to_prim_little_group = (
        []
    )  # [0..little_group_order) -> [0..prim_little_group_order)
    shifts = []  # (little_group_order, )
    for i in range(len(mapping_to_prim)):
        idx_prim = mapping_to_prim[i]  # in [0..order)
        idx_prim_little = remapping_prim_little_group.get(
            idx_prim
        )  # in [0..prim_little_group_order)
        if idx_prim_little is None:
            continue

        mapping_little_group.append(i)
        mapping_conv_to_prim_little_group.append(idx_prim_little)
        shifts.append(prim_translations[i] - uniq_prim_translations[idx_prim])

    # Adjust phase by centering translation
    phases = np.array([np.exp(-2j * np.pi * np.dot(prim_kpoint, shift)) for shift in shifts])
    irreps = []
    for prim_irrep in prim_irreps:
        # prim_irrep: (little_group_order, dim, dim)
        irrep = prim_irrep[mapping_conv_to_prim_little_group] * phases[:, None, None]
        irrep = purify_irrep_value(irrep)

        irreps.append(irrep)

    return irreps, rotations, translations, np.array(mapping_little_group)


def get_spacegroup_irreps_from_primitive_symmetry(
    rotations: NDArrayInt,
    translations: NDArrayFloat,
    kpoint: NDArrayFloat,
    real: bool = False,
    method: Literal["Neto", "random"] = "Neto",
    rtol: float = 1e-5,
    atol: float = 1e-8,
    max_num_random_generations: int = 4,
) -> tuple[list[NDArrayComplex], NDArrayInt] | tuple[list[NDArrayFloat], NDArrayInt]:
    r"""Compute all irreducible representations of given space group up to unitary transformation.

    Note that rotations and translations should be specified in a primitive cell.

    Parameters
    ----------
    rotations: array, (order, 3, 3)
        Assume a fractional coordinates `x` are transformed by the i-th symmetry operation as follows:
            ``np.dot(rotations[i, :, :], x) + translations[i, :]``
    translations: array, (order, 3)
    kpoint: array, (3, )
        Reciprocal vector with respect to reciprocal lattice.
        For pure translation :math:`\mathbf{t}`, returned irrep :math:`\Gamma^{(\alpha)}` takes

        .. math::
            \Gamma^{(\alpha)}((E, \mathbf{t})) = e^{ -i\mathbf{k}\cdot\mathbf{t} } \mathbf{1}.

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
        Relative tolerance
    atol: float
        Absolute tolerance to distinguish difference eigenvalues
    max_num_random_generations: int
        Maximum number of trials to generate random matrix

    Returns
    -------
    irreps: list of Irreps with (little_group_order, dim, dim)
        Let ``i = mapping_little_group[idx]``. ``irreps[alpha][i, :, :]`` is the ``alpha``-th irreducible matrix representation of ``(rotations[i], translations[i])``.
    mapping_little_group: array, (little_group_order, )
        Let ``i = mapping_little_group[idx]``.
        ``(rotations[i], translations[i])`` belongs to the little group of given space space group and kpoint.
    """
    # Sanity check to use primitive cell
    for rotation, translation in zip(rotations, translations):
        if np.allclose(rotation, np.eye(3), rtol=rtol, atol=atol) and not np.allclose(
            translation, 0, atol=atol
        ):
            raise ValueError("Specify symmetry operations in primitive cell!")

    little_rotations, little_translations, mapping_little_group = get_little_group(
        rotations, translations, kpoint, atol=atol
    )

    # Small representations of little group
    irreps, indicators = enumerate_small_representations(
        little_rotations,
        little_translations,
        kpoint,
        real=real,
        method=method,
        rtol=rtol,
        atol=atol,
        max_num_random_generations=max_num_random_generations,
    )

    return irreps, mapping_little_group


def get_crystallographic_pointgroup_irreps_from_symmetry(
    rotations: NDArrayInt,
    real: bool = False,
    method: Literal["Neto", "random"] = "Neto",
    rtol: float = 1e-5,
    atol: float = 1e-8,
    max_num_random_generations: int = 4,
) -> list[NDArrayComplex]:
    """Compute all irreducible representations of given crystallographic point group up to unitary transformation.

    Assume matrix representation of given crystallographic point group is in "standard" setting shown in Table 3.2.3.3 of International Table for Crystallography Vol. A (2016).

    Parameters
    ----------
    rotations: array, (order, 3, 3)
        Assume a point coordinates ``x`` are transformed into ``np.dot(rotations[i, :, :], x)`` by the ``i``-th symmetry operation.
    real: bool, default=False
        If True, return irreps over real vector space (so called physically irreducible representations)
    method: str, 'Neto' or 'random'
        'Neto': construct irreps from a fixed chain of subgroups of little co-group
        'random': construct irreps by numerically diagonalizing a random matrix commute with regular representation
    rtol: float
        Relative tolerance to distinguish difference eigenvalues
    atol: float
        Absolute tolerance to distinguish difference eigenvalues
    max_num_random_generations: int
        Maximum number of trials to generate random matrix

    Returns
    -------
    irreps: list of Irreps with (order, dim, dim)
    """
    irreps, indicators = enumerate_unitary_irreps(
        rotations,
        factor_system=None,
        real=real,
        method=method,
        rtol=rtol,
        atol=atol,
        max_num_random_generations=max_num_random_generations,
    )
    return irreps
