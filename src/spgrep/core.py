"""Enumerate unitary irreps from crystal structure or symmetry operations."""
from __future__ import annotations

from typing import Literal

import numpy as np
from spglib import get_magnetic_symmetry_dataset, get_symmetry_dataset

from spgrep.corep import enumerate_spinor_small_corepresentations
from spgrep.group import get_little_group
from spgrep.irreps import (
    enumerate_small_representations,
    enumerate_unitary_irreps,
    purify_irrep_value,
)
from spgrep.spinor import (
    enumerate_spinor_small_representations,
    get_spinor_factor_system,
)
from spgrep.transform import (
    get_primitive_transformation_matrix,
    transform_symmetry_and_kpoint,
    unique_primitive_symmetry,
)
from spgrep.utils import NDArrayBool, NDArrayComplex, NDArrayFloat, NDArrayInt

################################################################################
# Linear representation
################################################################################


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
    kpoint: array, (3, )
        Reciprocal vector with respect to ``reciprocal_lattice``
        For pure translation :math:`\mathbf{t}`, returned irrep :math:`\Gamma^{(\alpha)}` takes

        .. math::
            \Gamma^{(\alpha)}((E, \mathbf{t})) = e^{ -i\mathbf{k}\cdot\mathbf{t} } \mathbf{1}.

    method: str, 'Neto' or 'random'
        'Neto': construct irreps from a fixed chain of subgroups of little co-group
        'random': construct irreps by numerically diagonalizing a random matrix commute with regular representation
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
    rotations: array[int], (num_sym, 3, 3)
        Linear parts of symmetry operations
    translations: array, (num_sym, 3)
        Translation parts of symmetry operations
    mapping_little_group: array, (little_group_order, )
        Let ``i = mapping_little_group[idx]``.
        ``(rotations[i], translations[i])`` belongs to the little group of given space space group and kpoint.
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

    # Go back to conventional cell
    irreps, mapping_little_group = _adjust_phase_for_centering_translations(
        prim_translations,
        prim_kpoint,
        uniq_prim_translations,
        mapping_to_prim,
        prim_irreps,
        mapping_prim_little_group,
    )

    return irreps, rotations, translations, mapping_little_group


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

    Note that ``rotations`` and ``translations`` should be specified in a primitive cell.

    Parameters
    ----------
    rotations: array[int], (order, 3, 3)
        Assume a fractional coordinates `x` are transformed by the i-th symmetry operation as follows:
            ``np.dot(rotations[i, :, :], x) + translations[i, :]``
    translations: array, (order, 3)
    kpoint: array, (3, )
        Reciprocal vector with respect to reciprocal lattice.
        For pure translation :math:`\mathbf{t}`, returned irrep :math:`\Gamma^{(\alpha)}` takes

        .. math::
            \Gamma^{(\alpha)}((E, \mathbf{t})) = e^{ -i\mathbf{k}\cdot\mathbf{t} } \mathbf{1}.

        See :ref:`physically_irreps` for details.

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
    rotations: array[int], (order, 3, 3)
        Assume a point coordinates ``x`` are transformed into ``np.dot(rotations[i, :, :], x)`` by the ``i``-th symmetry operation.
    real: bool, default=False
        If True, return irreps over real vector space (so called physically irreducible representations).
        See :ref:`physically_irreps` for details.
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
    irreps, _ = enumerate_unitary_irreps(
        rotations,
        factor_system=None,
        real=real,
        method=method,
        rtol=rtol,
        atol=atol,
        max_num_random_generations=max_num_random_generations,
    )
    return irreps


################################################################################
# Spin representation
################################################################################


def get_spacegroup_spinor_irreps(
    lattice: NDArrayFloat,
    positions: NDArrayFloat,
    numbers: NDArrayInt,
    magmoms: NDArrayFloat | None = None,
    kpoint: NDArrayFloat | None = None,
    method: Literal["Neto", "random"] = "Neto",
    reciprocal_lattice: NDArrayFloat | None = None,
    symprec: float = 1e-5,
    rtol: float = 1e-5,
    atol: float = 1e-8,
    max_num_random_generations: int = 4,
) -> tuple[
    list[NDArrayComplex],
    NDArrayComplex,
    NDArrayComplex,
    NDArrayInt,
    NDArrayFloat,
    NDArrayInt,
] | tuple[
    list[NDArrayComplex],
    NDArrayComplex,
    NDArrayComplex,
    NDArrayBool,
    NDArrayInt,
    NDArrayFloat,
    NDArrayInt,
    NDArrayInt,
]:
    r"""Compute all irreducible representations :math:`\mathbf{\Gamma}^{\mathbf{k}\alpha}` of space group of given structure up to unitary transformation for spinor.

    Each irrep :math:`\mathbf{\Gamma}^{\mathbf{k}\alpha}` satisfies

    .. math::
       \mathbf{\Gamma}^{\mathbf{k}\alpha}((\mathbf{S}_{i}, \mathbf{w}_{i})) \mathbf{\Gamma}^{\mathbf{k}\alpha}((\mathbf{S}_{j}, \mathbf{w}_{j}))
       = z(\mathbf{S}_{i}, \mathbf{S}_{j}) \mathbf{\Gamma}^{\mathbf{k}\alpha}((\mathbf{S}_{i}, \mathbf{w}_{i})(\mathbf{S}_{j}, \mathbf{w}_{j})).

    See :ref:`spinor_factor_system` for Spgrep's convention of spin-derived factor system :math:`z(\mathbf{S}_{i}, \mathbf{S}_{j})`.

    Parameters
    ----------
    lattice: array, (3, 3)
        Row-wise basis vectors. ``lattice[i, :]`` is the i-th lattice vector.
    positions: array, (num_atoms, 3)
        Fractional coordinates of sites
    numbers: array, (num_atoms, )
        Integer list specifying atomic species
    magmoms: (Optional) array, (num_atoms, )
        Collinear magnetic moments. If specified, return co-representations.
        See :ref:{corep} for details.
    kpoint: array, (3, )
        Reciprocal vector with respect to ``reciprocal_lattice``
        For pure translation :math:`\mathbf{t}`, returned irrep :math:`\Gamma^{(\alpha)}` takes

        .. math::
            \Gamma^{(\alpha)}((E, \mathbf{t})) = e^{ -i\mathbf{k}\cdot\mathbf{t} } \mathbf{1}.

    method: str, 'Neto' or 'random'
        'Neto': construct irreps from a fixed chain of subgroups of little co-group
        'random': construct irreps by numerically diagonalizing a random matrix commute with regular representation
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
    little_spinor_factor_system: array, (little_group_order, little_group_order)
        ``spinor_factor_system[i, j]`` stands for factor system :math:`z(\mathbf{S}_{i}, \mathbf{S}_{j})`
    little_unitary_rotations: array, (little_group_order, 2, 2)
        SU(2) rotations on spinor.
    anti_linear: (Optional) array[bool], (little_group_order, )
        Appeared when ``time_reversal`` is specified.
        If ``anti_linear[i] == True``, the ``i``-th operator is anti-linear.
    rotations: array[int], (num_sym, 3, 3)
    translations: array, (num_sym, 3)
    time_reversals: array[int], (num_sym, )
    mapping_little_group: array, (little_group_order, )
        Let ``i = mapping_little_group[idx]``.
        (rotations[i], translations[i]) belongs to the little group of given space space group and kpoint.
    """
    if kpoint is None:
        kpoint = np.zeros(3)

    # Transform given `kpoint` in dual of `lattice`
    dual_lattice = np.linalg.inv(lattice).T
    if reciprocal_lattice is None:
        reciprocal_lattice = dual_lattice
    # kpoint @ reciprocal_lattice == kpoint_conv @ dual_lattice
    kpoint_conv = kpoint @ reciprocal_lattice @ np.linalg.inv(dual_lattice)

    if magmoms is None:
        dataset = get_symmetry_dataset(cell=(lattice, positions, numbers), symprec=symprec)
    else:
        dataset = get_magnetic_symmetry_dataset(
            cell=(lattice, positions, numbers, magmoms), symprec=symprec
        )

    rotations = dataset["rotations"]
    translations = dataset["translations"]

    # Transform to primitive
    to_primitive = get_primitive_transformation_matrix(dataset["hall_number"])
    prim_rotations, prim_translations, prim_kpoint = transform_symmetry_and_kpoint(
        to_primitive, rotations, translations, kpoint_conv
    )
    prim_lattice = to_primitive.T @ lattice  # (AP)^T = P^T @ A^T
    # mapping_to_prim: [0..num_sym) -> [0..order)
    uniq_prim_rotations, uniq_prim_translations, mapping_to_prim = unique_primitive_symmetry(
        prim_rotations, prim_translations
    )

    if magmoms is None:
        # mapping_prim_little_group: [0..prim_little_group_order) -> [0..order)
        (
            prim_irreps,
            little_spinor_factor_system,
            little_unitary_rotations,
            mapping_prim_little_group,
        ) = get_spacegroup_spinor_irreps_from_primitive_symmetry(  # type: ignore
            lattice=prim_lattice,
            rotations=uniq_prim_rotations,
            translations=uniq_prim_translations,
            kpoint=prim_kpoint,
            method=method,
            rtol=rtol,
            atol=atol,
            max_num_random_generations=max_num_random_generations,
        )
    else:
        time_reversals = dataset["time_reversals"]
        uniq_prim_time_reversals = time_reversals[mapping_to_prim]
        (
            prim_irreps,
            _,
            little_spinor_factor_system,
            little_unitary_rotations,
            little_anti_linear,
            mapping_prim_little_group,
        ) = get_spacegroup_spinor_irreps_from_primitive_symmetry(  # type: ignore
            lattice=prim_lattice,
            rotations=uniq_prim_rotations,
            translations=uniq_prim_translations,
            time_reversals=uniq_prim_time_reversals,
            kpoint=prim_kpoint,
            method=method,
            rtol=rtol,
            atol=atol,
            max_num_random_generations=max_num_random_generations,
        )

    # Go back to conventional cell
    irreps, mapping_little_group = _adjust_phase_for_centering_translations(
        prim_translations,
        prim_kpoint,
        uniq_prim_translations,
        mapping_to_prim,
        prim_irreps,
        mapping_prim_little_group,
    )

    if magmoms is None:
        return (
            irreps,
            little_spinor_factor_system,
            little_unitary_rotations,
            rotations,
            translations,
            mapping_little_group,
        )
    else:
        return (
            irreps,
            little_spinor_factor_system,
            little_unitary_rotations,
            little_anti_linear,
            rotations,
            translations,
            time_reversals,
            mapping_little_group,
        )


def get_spacegroup_spinor_irreps_from_primitive_symmetry(
    lattice: NDArrayFloat,
    rotations: NDArrayInt,
    translations: NDArrayFloat,
    time_reversals: NDArrayInt | None = None,
    kpoint: NDArrayFloat | None = None,
    method: Literal["Neto", "random"] = "Neto",
    rtol: float = 1e-5,
    atol: float = 1e-8,
    max_num_random_generations: int = 4,
) -> tuple[list[NDArrayComplex], NDArrayComplex, NDArrayComplex, NDArrayInt] | tuple[
    list[NDArrayComplex], list[int], NDArrayComplex, NDArrayComplex, NDArrayBool, NDArrayInt
]:
    r"""Compute all irreducible representations :math:`\mathbf{\Gamma}^{\mathbf{k}\alpha}` of given space group up to unitary transformation for spinor.

    Each irrep :math:`\mathbf{\Gamma}^{\mathbf{k}\alpha}` satisfies

    .. math::
       \mathbf{\Gamma}^{\mathbf{k}\alpha}((\mathbf{S}_{i}, \mathbf{w}_{i})) \mathbf{\Gamma}^{\mathbf{k}\alpha}((\mathbf{S}_{j}, \mathbf{w}_{j}))
       = z(\mathbf{S}_{i}, \mathbf{S}_{j}) \mathbf{\Gamma}^{\mathbf{k}\alpha}((\mathbf{S}_{i}, \mathbf{w}_{i})(\mathbf{S}_{j}, \mathbf{w}_{j})).

    Note that rotations and translations should be specified in a primitive cell.
    See :ref:`spinor_factor_system` for Spgrep's convention of spin-derived factor system :math:`z(\mathbf{S}_{i}, \mathbf{S}_{j})`.

    Parameters
    ----------
    lattice: array, (3, 3)
        Row-wise basis vectors. ``lattice[i, :]`` is the i-th lattice vector.
    rotations: array[int], (order, 3, 3)
        Assume a fractional coordinates `x` are transformed by the i-th symmetry operation as follows:
            ``np.dot(rotations[i, :, :], x) + translations[i, :]``
    translations: array, (order, 3)
    time_reversals: array[int] | None, (order, )
        If specified, return co-representations
        See :ref:{corep} for details.
    kpoint: array, (3, )
        Reciprocal vector with respect to reciprocal lattice.
        For pure translation :math:`\mathbf{t}`, returned irrep :math:`\Gamma^{(\alpha)}` takes

        .. math::
            \mathbf{\Gamma}^{\mathbf{k}\alpha}((E, \mathbf{t})) = e^{ -i\mathbf{k}\cdot\mathbf{t} } \mathbf{1}.
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
    indicators: (Optional) list[int]
        Appeared when ``time_reversal`` is specified.
        Frobenius-Schur indicators for each small representation
    little_spinor_factor_system: array, (little_group_order, little_group_order)
        ``spinor_factor_system[i, j]`` stands for factor system :math:`z(\mathbf{S}_{i}, \mathbf{S}_{j})`
    little_unitary_rotations: array, (little_group_order, 2, 2)
        SU(2) rotations on spinor.
    anti_linear: (Optional) array[bool], (order, )
        Appeared when ``time_reversal`` is specified.
        If ``anti_linear[i] == True``, the ``i``-th operator is anti-linear.
    mapping_little_group: array, (little_group_order, )
        Let ``i = mapping_little_group[idx]``.
        ``(rotations[i], translations[i])`` belongs to the little group of given space space group and kpoint.
    """
    if kpoint is None:
        kpoint = np.zeros(3)

    # Sanity check to use primitive cell
    for rotation, translation in zip(rotations, translations):
        if np.allclose(rotation, np.eye(3), rtol=rtol, atol=atol) and not np.allclose(
            translation, 0, atol=atol
        ):
            raise ValueError("Specify symmetry operations in primitive cell!")

    little_rotations, little_translations, mapping_little_group = get_little_group(
        rotations, translations, kpoint, atol=atol
    )

    if time_reversals is None:
        (
            irreps,
            little_spin_factor_system,
            little_unitary_rotations,
        ) = enumerate_spinor_small_representations(
            lattice=lattice,
            little_rotations=little_rotations,
            little_translations=little_translations,
            kpoint=kpoint,
            method=method,
            rtol=rtol,
            atol=atol,
            max_num_random_generations=max_num_random_generations,
        )

        return irreps, little_spin_factor_system, little_unitary_rotations, mapping_little_group
    else:
        # Co-representation
        little_time_reversals = time_reversals[mapping_little_group]
        (
            irreps,
            indicators,
            little_spin_factor_system,
            little_unitary_rotations,
            little_anti_linear,
        ) = enumerate_spinor_small_corepresentations(
            lattice,
            little_rotations,
            little_translations,
            little_time_reversals,
            kpoint=kpoint,
            method=method,
            rtol=rtol,
            atol=atol,
            max_num_random_generations=max_num_random_generations,
        )
        return (
            irreps,
            indicators,
            little_spin_factor_system,
            little_unitary_rotations,
            little_anti_linear,
            mapping_little_group,
        )


def get_crystallographic_pointgroup_spinor_irreps_from_symmetry(
    lattice: NDArrayFloat,
    rotations: NDArrayInt,
    time_reversals: NDArrayInt | None = None,
    method: Literal["Neto", "random"] = "Neto",
    rtol: float = 1e-5,
    atol: float = 1e-8,
    max_num_random_generations: int = 4,
) -> tuple[list[NDArrayComplex], NDArrayComplex, NDArrayComplex] | tuple[
    list[NDArrayComplex], list[int], NDArrayComplex, NDArrayComplex, NDArrayBool
]:
    r"""Compute all irreducible representations :math:`\mathbf{D}^{\alpha}` of given crystallographic point group up to unitary transformation for spinor.

    Each irrep :math:`\mathbf{D}^{\alpha}` satisfies

    .. math::
       \mathbf{D}^{\alpha}(\mathbf{S}_{i}) \mathbf{D}^{\alpha}(\mathbf{S}_{j})
       = z(\mathbf{S}_{i}, \mathbf{S}_{j}) \mathbf{D}^{\alpha}(\mathbf{S}_{k}).

    Assume matrix representation of given crystallographic point group is in "standard" setting shown in Table 3.2.3.3 of International Table for Crystallography Vol. A (2016).

    See :ref:`spinor_factor_system` for Spgrep's convention of spinor-derived factor system :math:`z(\mathbf{S}_{i}, \mathbf{S}_{j})`.

    Parameters
    ----------
    lattice: array, (3, 3)
        Row-wise basis vectors. ``lattice[i, :]`` is the i-th lattice vector.
    rotations: array[int], (order, 3, 3)
        Assume a point coordinates ``x`` are transformed into ``np.dot(rotations[i, :, :], x)`` by the ``i``-th symmetry operation.
    time_reversals: array[int] | None, (order, )
        If specified, return co-representations
        See :ref:{corep} for details.
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
    irreps: list of unitary irreps with (order, dim, dim)
    indicators: (Optional) list[int]
        Appeared when ``time_reversal`` is specified.
        Frobenius-Schur indicators for each small representation
    factor_system: array, (order, order)
        ``factor_system[i, j]`` stands for factor system :math:`z(\mathbf{S}_{i}, \mathbf{S}_{j})`
    unitary_rotations: array, (order, 2, 2)
        SU(2) rotations on spinor.
    anti_linear: (Optional) array[bool], (order, )
        Appeared when ``time_reversal`` is specified.
        If ``anti_linear[i] == True``, the ``i``-th operator is anti-linear.
    """
    if time_reversals is None:
        factor_system, unitary_rotations = get_spinor_factor_system(lattice, rotations)

        irreps, _ = enumerate_unitary_irreps(
            rotations,
            factor_system=factor_system,
            real=False,  # Nonsense to consider real-value irreps
            method=method,
            rtol=rtol,
            atol=atol,
            max_num_random_generations=max_num_random_generations,
        )
        return irreps, factor_system, unitary_rotations
    else:
        # Co-representation
        order = len(rotations)
        (
            irreps,
            indicators,
            factor_system,
            unitary_rotations,
            anti_linear,
        ) = enumerate_spinor_small_corepresentations(
            lattice,
            little_rotations=rotations,
            little_translations=np.zeros((order, 3)),
            little_time_reversals=time_reversals,
            kpoint=np.zeros(3),
            method=method,
            rtol=rtol,
            atol=atol,
            max_num_random_generations=max_num_random_generations,
        )
        return irreps, indicators, factor_system, unitary_rotations, anti_linear


################################################################################
# Auxiliary functions
################################################################################


def _adjust_phase_for_centering_translations(
    prim_translations,
    prim_kpoint,
    uniq_prim_translations,
    mapping_to_prim,
    prim_irreps,
    mapping_prim_little_group,
):
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

    phases = np.array([np.exp(-2j * np.pi * np.dot(prim_kpoint, shift)) for shift in shifts])

    irreps = []
    for prim_irrep in prim_irreps:
        # prim_irrep: (little_group_order, dim, dim)
        irrep = prim_irrep[mapping_conv_to_prim_little_group] * phases[:, None, None]
        irrep = purify_irrep_value(irrep)

        irreps.append(irrep)
    return irreps, np.array(mapping_little_group)
