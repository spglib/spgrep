"""Utility functions."""
from __future__ import annotations

from typing import Any, Literal

import numpy as np
from numpy.typing import NDArray
from spglib import get_symmetry_from_database
from typing_extensions import TypeAlias  # for Python<3.10

NDArrayBool: TypeAlias = NDArray[np.bool_]
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


def grassmann_distance(
    basis1: NDArrayComplex,
    basis2: NDArrayComplex,
    ord: Literal["min", "projection"] = "min",
) -> float:
    r"""Return Grassmann distance between two linear subspaces spanned by ``basis1`` and ``basis2``.

    References
    * [1] Jihun Hamm and Daniel D. Lee. 2008. Grassmann discriminant analysis: a unifying view on subspace-based learning. In Proceedings of the 25th international conference on Machine learning (ICML '08). Association for Computing Machinery, New York, NY, USA, 376–383. https://doi.org/10.1145/1390156.1390204
    * [2] Schubert Varieties and Distances between Subspaces of Different Dimensions, Ke Ye and Lek-Heng Lim, SIAM Journal on Matrix Analysis and Applications 2016 37:3, 1176-1197

    Parameters
    ----------
    basis1: array, (k, n)
        ``k``-dimensional subspace in :math:`\mathbb{C}^{n}`.
        ``basis1[i]`` is the ``i``-th basis vector.
    basis2: array, (l, n)
        ``l``-dimensional subspace in :math:`\mathbb{C}^{n}`.
        ``basis2[i]`` is the ``i``-th basis vector.
    ord: str
        Kind of Grassmann distance to be calculated
        * ``ord='min'``: Min correlation
        * ``ord='projection'``: Projection metric

    Returns
    -------
    distance: float
    """
    # Orthonormal bases
    # QR decomposition of column-wise vectors gives Gram-Schmidt orthonormalized vectors in column wise.
    col_orthonormal_basis1 = np.linalg.qr(np.transpose(basis1))[0]
    col_orthonormal_basis2 = np.linalg.qr(np.transpose(basis2))[0]

    # Singular values in descending order
    canonical_correlations = np.linalg.svd(
        np.conj(col_orthonormal_basis1.T) @ col_orthonormal_basis2, compute_uv=False
    )

    dim = min(basis1.shape[0], basis2.shape[0])
    if ord == "min":
        distance = np.sqrt(np.clip(1 - canonical_correlations[dim - 1] ** 2, a_min=0, a_max=None))
    elif ord == "projection":
        distance = np.sqrt(
            np.mean(np.arccos(np.clip(canonical_correlations[:dim], a_min=-1, a_max=1)) ** 2)
        )
    else:
        raise ValueError(f"Unknown type of Grassmann distance: {ord}")

    return distance


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
