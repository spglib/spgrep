from __future__ import annotations

from itertools import product
from typing import Any

import numpy as np
from numpy.typing import NDArray
from typing_extensions import TypeAlias  # for Python<3.10

NDArrayInt: TypeAlias = NDArray[np.int_]
NDArrayFloat: TypeAlias = NDArray[np.float_]


def is_matrix_group(rotations: NDArrayInt) -> bool:
    """
    Return True iff given integer matrices forms group.
    """
    rotations_set = {ndarray2d_to_integer_tuple(r) for r in rotations}

    for g1, g2 in product(rotations, repeat=2):
        # g1^-1 * g2 should be included in rotations_set
        g1_inv = np.linalg.inv(g1)
        if not is_integer_array(g1_inv):
            return False
        g1_inv = np.around(g1_inv).astype(int)
        g1_inv_g2 = np.dot(g1_inv, g2)
        if ndarray2d_to_integer_tuple(g1_inv_g2) not in rotations_set:
            return False

    return True


def is_integer_array(array: NDArrayFloat, rtol: float = 1e-5, atol: float = 1e-8) -> bool:
    array_int = np.around(array).astype(int)
    return np.allclose(array_int, array, rtol=rtol, atol=atol)


def ndarray2d_to_integer_tuple(array: NDArrayFloat) -> tuple[tuple[Any]]:
    array_int = np.around(array).astype(int)
    array_t = tuple(map(tuple, array_int.tolist()))
    return array_t  # type: ignore
