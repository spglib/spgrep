"""Utility functions."""
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
    """Return true if all values of ``array`` are almost integers."""
    array_int = np.around(array).astype(int)
    return np.allclose(array_int, array, rtol=rtol, atol=atol)


def ndarray2d_to_integer_tuple(array: NDArrayFloat) -> tuple[tuple[Any]]:
    """Convert two-dimensional array to tuple of tuple."""
    array_int = np.around(array).astype(int)
    array_t = tuple(map(tuple, array_int.tolist()))
    return array_t  # type: ignore


def get_symmetry_from_hall_number(hall_number: int) -> tuple[NDArrayInt, NDArrayFloat]:
    """Return symmetry operations from Hall number.

    Parameters
    ----------
    hall_number: int

    Returns
    -------
    rotations: (num_syms, 3, 3)
    translations: (num_syms, 3)
    """
    symmetry = get_symmetry_from_database(hall_number)
    rotations = symmetry["rotations"]
    translations = symmetry["translations"]
    return rotations, translations


def is_prime(n: int) -> bool:
    """Return true if given number is prime."""
    for i in range(2, n):
        if n % i == 0:
            return False
    return True


def nroot(z: np.complex_, n: int) -> np.complex_:
    """Return `n`-th power root of `z` with the minimum angle."""
    root = z ** (1 / n)
    r = np.absolute(root)
    angle = np.angle(root)
    angle -= np.rint(angle * n / (2 * np.pi)) * 2 * np.pi / n
    return r * np.exp(1j * angle)


def contain_space(
    basis1: NDArrayComplex,
    basis2: NDArrayComplex,
    atol: float = 1e-8,
) -> bool:
    """Return true if vector space spanned by ``basis2`` is contained in that by ``basis1``.

    That is, return True if any linear combination A[i, j] exists such that
        basis2[j] == sum_{i} basis1[i] * A[i, j]
    which is equivalent to ``A.T @ basis1 == basis2``.

    Parameters
    ----------
    basis1: array, (dim_irrep, dim)
    basis2: array, (dim_irrep, dim)
    """
    if len(basis1) == 0:
        return False  # basis1 is empty

    # Solve basis1.T @ A = basis2.T
    A, residual, _, _ = np.linalg.lstsq(basis1.T, basis2.T, rcond=None)

    # Always compare vectors by L_infinity norm
    basis2_T_near = basis1.T @ A
    return np.allclose(basis2_T_near, basis2.T, atol=atol)


def mode_dot(coeffs: NDArray, list_matrix: list[NDArray]) -> NDArray:
    """Calculate p-mode product of tensor and list of matrices.

    For example, 3-mode product of ``coeffs`` and ``list_matrix=[m1, m2, m3]`` is ``np.einsum("ijk,ia,jb,kc", coeffs, m1, m2, m3)``.

    Parameters
    ----------
    coeff: array, (m, ..., m)
    list_matrix: list of array, shape of the i-th array is (m, n_i)

    Returns
    -------
    ret: (n_1, ..., n_p)
    """

    def _mode_dot(tensor, mat, axis):
        # tensor: (m_{0}, ..., m_{axis}, ...)
        # mat: (m_{axis}, n)
        # ret: (m_{0}, ..., n, ...)
        dim_tensor = tuple(
            [
                1,
            ]
            + list(tensor.shape)
        )
        tensor_reshaped = tensor.reshape(dim_tensor)

        dim_mat = tuple(
            list(mat.shape)
            + [
                1,
            ]
            * (tensor.ndim - 1)
        )
        mat_reshaped = np.swapaxes(np.swapaxes(mat.reshape(dim_mat), 0, 1), 1, axis + 1)

        ret = np.swapaxes(
            np.sum(tensor_reshaped * mat_reshaped, axis=axis + 1, keepdims=True), 0, axis + 1
        )
        return ret[0]

    ret = coeffs.copy()
    for axis in range(len(list_matrix)):
        ret = _mode_dot(ret, list_matrix[axis], axis)
    return ret
