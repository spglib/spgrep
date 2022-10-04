"""Symmetric tensor."""
from __future__ import annotations

from itertools import permutations, product

import numpy as np

from spgrep.irreps import enumerate_unitary_irreps, is_equivalent_irrep
from spgrep.representation import get_character, get_direct_product, project_to_irrep
from spgrep.utils import NDArrayComplex, NDArrayFloat, NDArrayInt, contain_space


def get_symmetry_adapted_tensors(
    rep: NDArrayComplex | NDArrayFloat,
    rotations: NDArrayInt,
    rank: int,
    real: bool = False,
    atol: float = 1e-8,
) -> list[NDArrayComplex] | list[NDArrayFloat]:
    """Calculate symmetry-adapted tensors with rank=``rank``.

    Parameters
    ----------
    rep: array, (order, dim, dim)
        Representation matrices
    rotations: array, (order, 3, 3)
    rank: int
        Rank of returned tensor
    real: bool

    Returns
    -------
    tensors: list of symmetry-adapted ``rank``-tensor with (dim, ..., dim)
    """
    order = rep.shape[0]

    irreps, _ = enumerate_unitary_irreps(rotations, real=real, atol=atol)
    all_vector_basis = []
    for irrep in irreps:
        list_basis = project_to_irrep(rep, irrep, atol=atol)
        for basis in list_basis:
            all_vector_basis.append((basis, irrep))

    # Iteratively construct irreps for rank-`tensor` including not invariant ones
    all_tensor_basis = all_vector_basis[:]
    for _ in range(rank - 1):
        next_all_tensor_basis = []
        for (tensor_basis, tensor_irrep), (vector_basis, vector_irrep) in product(
            all_tensor_basis, all_vector_basis
        ):
            # tensor_basis: (dim_tensor, dim(1), ..., dim(p-1))
            # vector_basis: (dim_vector, dim)
            dim_tensor = tensor_irrep.shape[1]
            dim_vector = vector_irrep.shape[1]
            direct_rep = get_direct_product(tensor_irrep, vector_irrep)
            for irrep in irreps:
                dim_irrep = irrep.shape[1]
                decomposed = project_to_irrep(direct_rep, irrep, atol=atol)
                for coeff in decomposed:
                    coeff = coeff.reshape(dim_irrep, dim_tensor, dim_vector)
                    # new_basis: (dim_irrep, dim(1), ..., dim(p))
                    new_basis = np.einsum(
                        "ijk,j...,ka->i...a", coeff, tensor_basis, vector_basis, optimize="greedy"
                    )
                    next_all_tensor_basis.append((new_basis, irrep))

        all_tensor_basis = next_all_tensor_basis

    # Take invariant tensors
    identity_character = np.ones((order,))
    tensors = []
    for basis, irrep in all_tensor_basis:
        if irrep.shape[1] > 1:
            continue
        character = get_character(irrep)
        if is_equivalent_irrep(character, identity_character):
            tensors.append(basis[0, ...])

    return tensors


def apply_intrinsic_symmetry(
    tensors: list[NDArrayComplex] | list[NDArrayFloat],
    atol: float = 1e-8,
) -> list[NDArrayComplex] | list[NDArrayFloat]:
    """Apply symmetric group on tensors.

    Note
    ----
    Current implementation may return wrong number of symmetrized tensors for higher rank...
    """
    list_sym_basis = []  # type: ignore
    for coeff in tensors:
        rank = len(coeff.shape)
        perms = list(permutations(range(rank)))
        # Apply intrinsic symmetry
        sym_tensor = np.zeros_like(coeff)
        coeff = coeff / np.linalg.norm(coeff)
        for perm in perms:
            sym_tensor += np.transpose(coeff, axes=perm)
        if np.allclose(sym_tensor, 0, atol=atol):
            continue
        sym_tensor /= np.linalg.norm(sym_tensor)

        # Store if independent to other symmetric tensors
        if len(list_sym_basis) > 0 and contain_space(
            np.concatenate(list_sym_basis, axis=0),
            sym_tensor.reshape(1, -1),
            atol=atol,
        ):
            continue

        list_sym_basis.append(sym_tensor.reshape(1, -1))

    # Reshape tensors
    sym_tensors = [basis.reshape(tensors[0].shape) for basis in list_sym_basis]

    return sym_tensors
