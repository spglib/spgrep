from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray
from spglib import get_symmetry_from_database
from typing_extensions import TypeAlias  # for Python<3.10

NDArrayInt: TypeAlias = NDArray[np.int_]
NDArrayFloat: TypeAlias = NDArray[np.float_]
NDArrayComplex: TypeAlias = NDArray[np.complex_]


def is_integer_array(array: NDArrayFloat, rtol: float = 1e-5, atol: float = 1e-8) -> bool:
    array_int = np.around(array).astype(int)
    return np.allclose(array_int, array, rtol=rtol, atol=atol)


def ndarray2d_to_integer_tuple(array: NDArrayFloat) -> tuple[tuple[Any]]:
    array_int = np.around(array).astype(int)
    array_t = tuple(map(tuple, array_int.tolist()))
    return array_t  # type: ignore


def get_symmetry_from_hall_number(hall_number: int) -> tuple[NDArrayInt, NDArrayFloat]:
    symmetry = get_symmetry_from_database(hall_number)
    rotations = symmetry["rotations"]
    translations = symmetry["translations"]
    return rotations, translations
