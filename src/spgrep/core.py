from __future__ import annotations

from spglib import get_symmetry

from spgrep.utils import NDArrayFloat, NDArrayInt


def get_spacegroup_irreps(
    lattice: NDArrayFloat,
    positions: NDArrayFloat,
    numbers: NDArrayInt,
    kpoint: NDArrayFloat,
    symprec: float = 1e-5,
):
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
        Reciprocal vector with respect to reciprocal lattice of `lattice`
    symprec: float

    Returns
    -------
    TODO
    """
    cell = (lattice, positions, numbers)
    symmetry = get_symmetry(cell, symprec=symprec)
    return get_spacegroup_irreps_from_symmetry(
        rotations=symmetry["rotations"],
        translations=symmetry["translations"],
        kpoint=kpoint,
    )


def get_spacegroup_irreps_from_symmetry(
    rotations: NDArrayInt, translations: NDArrayFloat, kpoint: NDArrayFloat
):
    """Compute all irreducible representations of given space group up to unitary transformation.

    Parameters
    ----------
    rotations: array, (num_sym, 3, 3)
        Assume a fractional coordinates `x` are transformed by the i-th symmetry operation as follows:
            np.dot(rotations[i, :, :], x) + translations[i, :]
        Note that rotations and translations can be specified in conventional cell.
    translations: array, (num_sym, 3)
    kpoint: array, (3, )
        Reciprocal vector with respect to reciprocal lattice of `lattice`

    Returns
    -------
    TODO
    """
    raise NotImplementedError


def get_crystallographic_pointgroup_irreps_from_symmetry(
    rotations: NDArrayInt,
):
    """Compute all irreducible representations of given crystallographic point group up to unitary transformation.
    Assume matrix representation of given crystallographic point group is in "standard" setting shown in Table 3.2.3.3 of International Table for Crystallography Vol. A (2016).

    Parameters
    ----------
    rotations: array, (order, 3, 3)
        Assume a point coordinates `x` are transformed into `np.dot(rotations[i, :, :], x)` by the i-th symmetry operation.

    Returns
    -------
    TODO
    """
    raise NotImplementedError
