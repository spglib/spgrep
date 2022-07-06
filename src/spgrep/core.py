from __future__ import annotations

import numpy as np
from spglib import get_symmetry_dataset

from spgrep.group import get_factor_system_from_little_group, get_little_group
from spgrep.irreps import (
    get_irreps,
    get_projective_regular_representation,
    get_regular_representation,
)
from spgrep.transform import (
    get_primitive_transformation_matrix,
    transform_symmetry_and_kpoint,
)
from spgrep.utils import NDArrayComplex, NDArrayFloat, NDArrayInt


def get_spacegroup_irreps(
    lattice: NDArrayFloat,
    positions: NDArrayFloat,
    numbers: NDArrayInt,
    kpoint: NDArrayFloat,
    reciprocal_lattice: NDArrayFloat | None = None,
    symprec: float = 1e-5,
) -> list[NDArrayComplex]:
    """Compute all irreducible representations of space group of given structure up to unitary transformation.

    Parameters
    ----------
    lattice: array, (3, 3)
        Row-wise basis vectors. lattice[i, :] is the i-th lattice vector.
    positions: array, (num_atoms, 3)
        Fractional coordinates of sites
    numbers: array, (num_atoms, )
        Integer list specifying atomic species
    kpoint: array, (3, )
        Reciprocal vector with respect to ``reciprocal_lattice``
    reciprocal_lattice: (Optional) array, (3, 3)
        ``reciprocal_lattice[i, :]`` is the i-th basis vector of reciprocal lattice for ``kpoint`` without `2 * pi factor`.
        If not specified, ``reciprocal_lattice`` is set to ``np.linalg.inv(lattice).T``.
    symprec: float

    Returns
    -------
    irreps: list of Irreps with (order, dim, dim)
    """
    # Transform given `kpoint` in dual of `lattice`
    dual_lattice = np.linalg.inv(lattice).T
    if reciprocal_lattice is None:
        reciprocal_lattice = dual_lattice
    # kpoint @ reciprocal_lattice == kpoint_conv @ dual_lattice
    kpoint_conv = kpoint @ reciprocal_lattice @ np.linalg.inv(dual_lattice)

    dataset = get_symmetry_dataset(cell=(lattice, positions, numbers), symprec=symprec)

    to_primitive = get_primitive_transformation_matrix(dataset["hall_number"])
    primitive_rotations, primitive_translations, primitive_kpoint = transform_symmetry_and_kpoint(
        to_primitive, dataset["rotations"], dataset["translations"], kpoint_conv
    )

    primitive_irreps = get_spacegroup_irreps_from_primitive_symmetry(
        rotations=primitive_rotations,
        translations=primitive_translations,
        kpoint=primitive_kpoint,
    )

    # TODO
    raise NotImplementedError
    return primitive_irreps


def get_spacegroup_irreps_from_primitive_symmetry(
    rotations: NDArrayInt,
    translations: NDArrayFloat,
    kpoint: NDArrayFloat,
    rtol: float = 1e-5,
    max_num_random_generations: int = 4,
) -> tuple[list[NDArrayComplex], NDArrayInt]:
    """Compute all irreducible representations of given space group up to unitary transformation.
    Note that rotations and translations should be specified in a primitive cell.

    Parameters
    ----------
    rotations: array, (order, 3, 3)
        Assume a fractional coordinates `x` are transformed by the i-th symmetry operation as follows:
            np.dot(rotations[i, :, :], x) + translations[i, :]
    translations: array, (order, 3)
    kpoint: array, (3, )
        Reciprocal vector with respect to reciprocal lattice
    rtol: float
        Relative tolerance
    max_num_random_generations: int
        Maximal number of trials to generate random matrix

    Returns
    -------
    irreps: list of Irreps with (little_group_order, dim, dim)
        Let ``i = mapping_little_group[idx]``. ``irreps[alpha][i, :, :]`` is the ``alpha``-th irreducible matrix representation of ``(rotations[i], translations[i])``.
    mapping_little_group: array, (little_group_order, )
        Let ``i = mapping_little_group[idx]``.
        (rotations[i], translations[i]) belongs to the little group of given space space group and kpoint.
    """
    # Sanity check to use primitive cell
    for rotation, translation in zip(rotations, translations):
        if np.allclose(rotation, np.eye(3), rtol=rtol) and not np.allclose(
            translation, 0, rtol=rtol
        ):
            raise ValueError("Specify symmetry operations in primitive cell!")

    little_rotations, little_translations, mapping_little_group = get_little_group(
        rotations, translations, kpoint, rtol
    )
    factor_system = get_factor_system_from_little_group(
        little_rotations, little_translations, kpoint
    )
    reg = get_projective_regular_representation(little_rotations, factor_system)
    small_reps = get_irreps(reg, rtol, max_num_random_generations)

    irreps = []
    for rep in small_reps:
        phases = np.array(
            [
                np.exp(-2j * np.pi * np.dot(kpoint, translation))
                for translation in little_translations
            ]
        )
        irreps.append(rep * phases[:, None, None])

    # TODO: symmetrize irreps
    return irreps, mapping_little_group


def get_crystallographic_pointgroup_irreps_from_symmetry(
    rotations: NDArrayInt,
    rtol: float = 1e-5,
    max_num_random_generations: int = 4,
) -> list[NDArrayComplex]:
    """Compute all irreducible representations of given crystallographic point group up to unitary transformation.
    Assume matrix representation of given crystallographic point group is in "standard" setting shown in Table 3.2.3.3 of International Table for Crystallography Vol. A (2016).

    Parameters
    ----------
    rotations: array, (order, 3, 3)
        Assume a point coordinates `x` are transformed into `np.dot(rotations[i, :, :], x)` by the i-th symmetry operation.
    rtol: float
        Relative tolerance to distinguish difference eigenvalues
    max_num_random_generations: int
        Maximal number of trials to generate random matrix

    Returns
    -------
    irreps: list of Irreps with (order, dim, dim)
    """
    reg = get_regular_representation(rotations)
    irreps = get_irreps(reg, rtol, max_num_random_generations)
    # TODO: symmetrize irreps
    return irreps
