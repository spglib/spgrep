from __future__ import annotations

import numpy as np

from spgrep._constants import ATOL, MAX_NUM_RANDOM_GENERATIONS
from spgrep.utils import NDArrayComplex, NDArrayFloat, nroot

from .representation import get_intertwiner


def get_physically_irrep(
    irrep: NDArrayComplex,
    indicator: int,
    atol: float = ATOL,
    max_num_random_generations: int = MAX_NUM_RANDOM_GENERATIONS,
) -> NDArrayFloat:
    """Compute physically irreducible representation (over real number) from given unitary irrep over complex number.

    Parameters
    ----------
    irrep: array, (order, dim, dim)
        Unitary (projective) irrep
    indicator: int
        Frobenius-Schur indicator: -1, 0, or 1
    atol: float
        Relative tolerance to compare
    max_num_random_generations: int
        Maximum number of trials to generate random matrix

    Returns
    -------
    real_irrep: array, (order, dim2, dim2)
        Physically irreducible representation composed of `irrep` and its conjugated irrep
        When `indicator==1`, `dim2 == dim`.
        When `indicator==-1` or `indicator==0`, `dim2 == 2 * dim`.
    """
    order = irrep.shape[0]
    dim = irrep.shape[1]

    # Assume kpoint is commensurated
    if indicator == 1:
        # Intertwiner with determinant=1
        conj_irrep = np.conj(irrep)
        U = get_intertwiner(
            irrep, conj_irrep, atol=atol, max_num_random_generations=max_num_random_generations
        )

        # Take real or imaginary part of eigenvectors for new basis vectors
        eigvals, eigvecs = np.linalg.eig(U)  # eigvecs[:, i] is the i-th eigenvector
        real_eigvecs = []
        for eigvec in np.transpose(eigvecs):
            if not np.allclose(np.real(eigvec), 0, atol=atol):
                real_eigvec = np.real(eigvec)
            else:
                real_eigvec = np.imag(eigvec)
            real_eigvecs.append(real_eigvec / np.linalg.norm(real_eigvec))
        S = np.array(real_eigvecs)

        # Square root of intertwiner
        T = S.T @ np.diag([nroot(eigval, 2) for eigval in eigvals]) @ S
        assert np.allclose(T @ T, U, atol=atol), "T is not square root of intertwiner."

        real_irrep = np.real(np.einsum("il,klm,mj->kij", np.conj(T), irrep, T, optimize="greedy"))

    elif indicator in [-1, 0]:
        real_irrep = np.empty((order, 2 * dim, 2 * dim), dtype=np.float64)
        # [ [Re D(g),  Im D(g)]
        #   [-Im D(g), Re D(g)] ]
        real_irrep[:, :dim, :dim] = np.real(irrep)
        real_irrep[:, :dim, dim:] = np.imag(irrep)
        real_irrep[:, dim:, :dim] = -np.imag(irrep)
        real_irrep[:, dim:, dim:] = np.real(irrep)

    return real_irrep
