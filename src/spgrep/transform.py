from typing import Literal

import numpy as np
from spglib import get_spacegroup_type

from spgrep.utils import NDArrayFloat, NDArrayInt


def transform_symmetry_and_kpoint(
    transformation_matrix: NDArrayFloat,
    rotations: NDArrayInt,
    translations: NDArrayFloat,
    kpoint: NDArrayFloat,
) -> tuple[NDArrayFloat, NDArrayFloat, NDArrayFloat]:
    """Return symmetry and k-vector coefficients in transformed coordinates.

    Let a given transformation_matrix be ``P``.
    Symmetry operation (R, t) is transformed to
        (P, 0)^-1 (R, t) (P, 0) = (P^-1 R P, P^-1 t).
    Coefficients of k-vector are transformed to P.T @ kpoint.
    """

    # (R, t) -> (P^-1 R P, P^-1 t)
    pinv = np.linalg.inv(transformation_matrix)
    transformed_rotations = np.array(
        [pinv @ rotation @ transformation_matrix for rotation in rotations]
    )
    transformed_translations = np.array([pinv @ translation for translation in translations])
    # k -> P^T k
    transformed_kpoint = transformation_matrix.T @ kpoint

    return transformed_rotations, transformed_translations, transformed_kpoint


def get_primitive_transformation_matrix(hall_number: int) -> NDArrayFloat:
    """Return transformation matrix from standard unit cell specified with hall_number into a primitive cell. The transformation matrix is consistent with Spglib's convention [1,2] and KVEC's convention [3].

    [1] https://spglib.github.io/spglib/definition.html
    [2] https://github.com/spglib/spglib/pull/137
    [3] M. I. Aroyo, D. Orobengoa, G. de la Flor, E.S. Tasci, J. M. Perez-Mato and H. Wondratschek, Acta Cryst. A70 126-137 (2014).
    """
    crystal_system = get_crystal_system(hall_number)
    spacegroup_type = get_spacegroup_type(hall_number)

    if crystal_system == "trigonal":
        choice = spacegroup_type["choice"]
        if choice == "H":
            return np.array(
                [
                    [2 / 3, -1 / 3, -1 / 3],
                    [1 / 3, 1 / 3, -2 / 3],
                    [1 / 3, 1 / 3, 1 / 3],
                ]
            )
        elif choice in ["", "R"]:
            return np.eye(3)
        else:
            raise ValueError("Unreachable!")

    hall_symbol = spacegroup_type["hall_symbol"]
    centering = get_centering(hall_symbol)
    if centering == "P":
        return np.eye(3)
    elif centering == "A":
        return np.array(
            [
                [1, 0, 0],
                [0, 1 / 2, -1 / 2],
                [0, 1 / 2, 1 / 2],
            ]
        )
    elif centering == "C":
        return np.array(
            [
                [1 / 2, 1 / 2, 0],
                [-1 / 2, 1 / 2, 0],
                [0, 0, 1],
            ]
        )
    elif centering == "I":
        return np.array(
            [
                [-1 / 2, 1 / 2, 1 / 2],
                [1 / 2, -1 / 2, 1 / 2],
                [1 / 2, 1 / 2, -1 / 2],
            ]
        )
    elif centering == "F":
        return np.array(
            [
                [0, 1 / 2, 1 / 2],
                [1 / 2, 0, 1 / 2],
                [1 / 2, 1 / 2, 0],
            ]
        )

    raise ValueError("Unreachable!")


def get_crystal_system(
    hall_number: int,
) -> Literal[
    "triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic"
]:
    crystal_system_range = {
        "triclinic": [1, 2],
        "monoclinic": [3, 107],
        "orthorhombic": [108, 348],
        "tetragonal": [349, 429],
        "trigonal": [430, 461],
        "hexagonal": [462, 488],
        "cubic": [489, 530],
    }

    for csystem, (lb, ub) in crystal_system_range.items():
        if (lb <= hall_number) and (hall_number <= ub):
            return csystem  # type: ignore

    raise ValueError("Unknown Hall number: {}".format(hall_number))


def get_centering(hall_symbol: str) -> Literal["P", "A", "C", "I", "R", "F"]:
    if hall_symbol[0] == "-":
        # e.g. "-P 2c 2b"
        return hall_symbol[1]  # type: ignore
    else:
        return hall_symbol[0]  # type: ignore
