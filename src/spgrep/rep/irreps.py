from __future__ import annotations

from itertools import product
from warnings import warn

import numpy as np

from spgrep.utils import (
    NDArrayComplex,
    grassmann_distance,
)

from .representation import get_character


def project_to_irrep(
    representation: NDArrayComplex,
    irrep: NDArrayComplex,
    atol: float = 1e-6,  # A bit large tolerance setting to handle numerical noise in `representation`
    max_num_trials: int = 10,
) -> list[NDArrayComplex]:
    r"""Construct basis functions for ``irrep`` by linear combinations of basis functions of ``representation``.

    Parameters
    ----------
    representation: array, (order, dim, dim)
    irrep: array, (order, dim_irrep, dim_irrep)
        Unitary (projective) irrep with factor system s.t. :math:`\mu(E, E) = 1`.
    atol: float, default=1e-5
        Absolute tolerance to compare basis vectors
    max_num_trials: int, default=10
        Maximum number to retry when failed to select projected basis vectors

    Returns
    -------
    basis: list of array with (irrep_dim, dim)
        Each basis vectors are orthonormal.
    """
    order = irrep.shape[0]
    dim_irrep = irrep.shape[1]
    dim = representation.shape[1]
    if representation.shape != (order, dim, dim) or irrep.shape != (order, dim_irrep, dim_irrep):
        raise ValueError("Given representation and irrep do not have consistent dimensions.")

    # Pre-compute number of independent basis vectors
    character_irrep = get_character(irrep)
    character = get_character(representation)
    character_sum = np.sum(np.conj(character_irrep) * character)
    if not np.isclose(character_sum, np.around(character_sum), atol=atol):
        warn("Inner product of characters should return an integer.")
    num_basis = np.around(character_sum) / order
    num_basis = np.around(np.real(num_basis)).astype(int)
    if num_basis == 0:
        return []

    def _project_to_irrep(adjusted_atol):
        count = 0
        basis: list[NDArrayComplex] = []
        for n in range(dim):
            for j in range(dim_irrep):
                # basis_nj[i, :] is the i-th basis vector forms given irrep (i = 0, ... dim_irrep-1)
                # These basis vectors are mutually orthogonal by construction!
                basis_nj = (
                    dim_irrep
                    / order
                    * np.einsum(
                        "ki,km->im",
                        np.conj(irrep[:, :, j]),
                        representation[:, :, n],
                        optimize="greedy",
                    )
                )

                norms = np.linalg.norm(basis_nj, axis=1)
                if np.any(np.isclose(norms, 0, atol=adjusted_atol)):
                    continue

                # Normalize basis vectors s.t. they are orthonormal.
                basis_nj /= norms[:, None]

                # Check if linearly independent with other basis vectors
                # If basis_nj is not independent, Grassmann distance (min correlation) should be one.
                # We use very rough tolerance, 0.5 to avoid numerical noises.
                if (len(basis) > 0) and grassmann_distance(
                    basis_nj, np.concatenate(np.array(basis), axis=0)
                ) < 0.5:
                    continue

                basis.append(basis_nj)
                count += 1

        return basis, count

    # Binary search for appropriate tolerance
    atol_lb = 1e-10
    atol_ub = 1e-2
    adjusted_atol = np.clip(atol, atol_lb, atol_ub)
    for _ in range(max_num_trials):
        basis, count = _project_to_irrep(adjusted_atol)
        if count == num_basis:
            break
        elif count < num_basis:
            # Tighten tolerance to compare basis vectors
            adjusted_atol = np.sqrt(atol_lb * adjusted_atol)
            warn(f"Tighten tolerance for projection: {adjusted_atol}")
        else:
            # Loosen tolerance to compare basis vectors
            adjusted_atol = np.sqrt(atol_ub * adjusted_atol)
            warn(f"Loosen tolerance for projection: {adjusted_atol}")

    if count < num_basis:
        warn(
            f"Inconsistent number of independent basis vectors (expect={num_basis}, actual={count})."
            "Try decreasing atol."
        )
    elif count > num_basis:
        warn(
            f"Inconsistent number of independent basis vectors (expect={num_basis}, actual={count})."
            "Try increasing atol."
        )

    return basis


def frobenius_schur_indicator(irrep: NDArrayComplex) -> int:
    r"""Inspect given unitary (projective) irrep is real, pseudo-real, or not unitary equivalent.

    .. math::
       \mathrm{indicator} =
       \frac{1}{|G|} \sum_{ g \in G } \chi(g^{2})

    Parameters
    ----------
    irrep: array, (order, dim, dim)

    Returns
    -------
    indicator: int
        If ``indicator==1``, it is real Reps.
        If ``indicator==-1``, it is psedu-real Reps.
        Otherwise, it and adjoint Reps. are not equivalent.
    """
    order = irrep.shape[0]
    indicator = np.einsum("kij,kji->", irrep, irrep, optimize="greedy") / order
    indicator = int(np.around(np.real(indicator)))

    if indicator > 1:
        raise ValueError(f"Given representation is not irreducible: indicator={indicator}")

    return indicator


def is_equivalent_irrep(character1: NDArrayComplex, character2: NDArrayComplex) -> bool:
    """Return true if two irreps are equivalent."""
    order = character1.shape[0]
    if np.around(np.sum(np.conj(character1) * character2)) == order:
        return True
    else:
        return False


def is_unique_irreps(irreps: list[NDArrayComplex]):
    characters = [get_character(irrep) for irrep in irreps]
    for (i, ci), (j, cj) in product(enumerate(characters), repeat=2):
        if is_equivalent_irrep(ci, cj) != (i == j):
            return False
    return True
