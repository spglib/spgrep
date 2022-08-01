"""Transformation of space group and kpoints."""
from __future__ import annotations

from typing import Literal

import numpy as np
from spglib import get_spacegroup_type

from spgrep.utils import NDArrayFloat, NDArrayInt, ndarray2d_to_integer_tuple


def transform_symmetry_and_kpoint(
    transformation_matrix: NDArrayFloat,
    rotations: NDArrayInt,
    translations: NDArrayFloat,
    kpoint: NDArrayFloat,
) -> tuple[NDArrayInt, NDArrayFloat, NDArrayFloat]:
    """Return symmetry and k-vector coefficients in transformed coordinates.

    This function does not unique duplicated symmetry operations after applying transformation_matrix.

    Let a given transformation_matrix be ``P``.
    Symmetry operation :math:`(R, t)` is transformed to
    :math:`(P, 0){^-1} (R, t) (P, 0) = (P^{-1} R P, P^{-1} t)`.
    Coefficients of k-vector are transformed to ``P.T @ kpoint``.

    Parameters
    ----------
    transformation_matrix: array, (3, 3)
    rotations: array, (num_sym, 3, 3)
    translations: array, (num_sym, 3)
    kpoint: array, (3, )

    Returns
    -------
    transformed_rotations: array, (num_sym, 3, 3)
    transformed_translations: array, (num_sym, 3)
    transformed_kpoint: array, (3, )
    """
    # (R, t) -> (P^-1 R P, P^-1 t)
    pinv = np.linalg.inv(transformation_matrix)
    transformed_rotations = np.array(
        [pinv @ rotation @ transformation_matrix for rotation in rotations]
    )
    transformed_rotations = np.around(transformed_rotations).astype(int)
    # Never take modulus for translations here!
    transformed_translations = np.array([pinv @ translation for translation in translations])
    # k -> P^T k
    transformed_kpoint = transformation_matrix.T @ kpoint

    return transformed_rotations, transformed_translations, transformed_kpoint


def unique_primitive_symmetry(
    rotations: NDArrayInt, translations: NDArrayFloat
) -> tuple[NDArrayInt, NDArrayFloat, list[int]]:
    """Remove duplicated symmetry operations.

    Parameters
    ----------
    rotations: array, (num_sym, 3, 3)
    translations: array, (num_sym, 3)

    Returns
    -------
    unique_rotations: array, (new_num_sym, 3, 3)
    unique_translations: array, (new_num_sym, 3)
    mapping_to_primitive_symmetry: array, (num_sym, )
        ``(rotations[i], translations[i])`` is transformed to ``(unique_rotations[j], unique_translations[j])`` where ``j = mapping_to_primitive_symmetry[i]``.
    """
    rotations_int: list[tuple[tuple[int]]] = []

    unique_rotations = []
    unique_translations = []
    mapping_to_primitive_symmetry = [-1 for _ in range(len(rotations))]
    for i, (rotation, translation) in enumerate(zip(rotations, translations)):
        rotation_int = ndarray2d_to_integer_tuple(rotation)
        try:
            j = rotations_int.index(rotation_int)
            # Duplicated symmetry operation
            mapping_to_primitive_symmetry[i] = j
        except ValueError:
            # New symmetry operation
            unique_rotations.append(rotation)
            unique_translations.append(np.remainder(translation, 1))
            mapping_to_primitive_symmetry[i] = i
            rotations_int.append(rotation_int)

    return np.array(unique_rotations), np.array(unique_translations), mapping_to_primitive_symmetry


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
    elif centering == "B":
        return np.array(
            [
                [1 / 2, 0, -1 / 2],
                [0, 1, 0],
                [1 / 2, 0, 1 / 2],
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
    """Return crystal system from Hall number."""
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
    """Return centering symbol from Hall symbol."""
    if hall_symbol[0] == "-":
        # e.g. "-P 2c 2b"
        return hall_symbol[1]  # type: ignore
    else:
        return hall_symbol[0]  # type: ignore
