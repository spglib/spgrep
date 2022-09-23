from __future__ import annotations

import numpy as np
import pytest
from spglib import get_symmetry_from_database

from spgrep.group import get_cayley_table
from spgrep.representation import is_representation
from spgrep.tensor import apply_intrinsic_symmetry, get_symmetry_adapted_tensors


def get_standard_basis() -> list[np.ndarray]:
    # Basis for symmetric matrix in Voigt order (xx, yy, zz, yz, zx, xy)
    basis = [
        np.array(
            [
                [1, 0, 0],
                [0, 0, 0],
                [0, 0, 0],
            ],
            dtype=np.float_,
        ),
        np.array(
            [
                [0, 0, 0],
                [0, 1, 0],
                [0, 0, 0],
            ],
            dtype=np.float_,
        ),
        np.array(
            [
                [0, 0, 0],
                [0, 0, 0],
                [0, 0, 1],
            ],
            dtype=np.float_,
        ),
        np.array(
            [
                [0, 0, 0],
                [0, 0, 1],
                [0, 1, 0],
            ],
            dtype=np.float_,
        )
        / np.sqrt(2),
        np.array(
            [
                [0, 0, 1],
                [0, 0, 0],
                [1, 0, 0],
            ],
            dtype=np.float_,
        )
        / np.sqrt(2),
        np.array(
            [
                [0, 1, 0],
                [1, 0, 0],
                [0, 0, 0],
            ],
            dtype=np.float_,
        )
        / np.sqrt(2),
    ]
    return basis


def get_representation_on_symmetric_matrix(rotations: np.ndarray) -> np.ndarray:
    # take [e_{1,1}, e_{2,2}, e_{3,3}, e_{2,3}, e_{3,1}, e_{1,2}] as basis
    basis = get_standard_basis()
    rep = np.zeros((len(rotations), len(basis), len(basis)))
    for pos, rotation in enumerate(rotations):
        for j, bj in enumerate(basis):
            # operated = rotation.T @ bj @ rotation
            operated = rotation @ bj @ rotation.T
            for i, bi in enumerate(basis):
                rep[pos, i, j] = np.sum(operated * bi) / np.sum(bi * bi)

    # Sanity check if `rep` satisfy property of representation
    table = get_cayley_table(rotations)
    assert is_representation(rep, table)

    return rep


@pytest.mark.parametrize(
    "hall_number,rank,num_expect",
    [
        # Pm-3m (No. 221)
        (517, 1, 1),
        (517, 2, 3),
        (517, 3, 6),
        (517, 4, 11),
        # P6/mmm (No. 191)
        (485, 1, 2),
        (485, 2, 5),
        (485, 3, 10),
        # (485, 4, 18),  # TODO
    ],
)
def test_symmetric_tensor(hall_number, rank, num_expect):
    symmetry = get_symmetry_from_database(hall_number=hall_number)
    rotations = symmetry["rotations"]

    rep = get_representation_on_symmetric_matrix(rotations)
    atol = 1e-8
    tensors = get_symmetry_adapted_tensors(rep, rotations, rank, real=True, atol=atol)
    sym_tensors = apply_intrinsic_symmetry(tensors, atol=atol)

    assert len(sym_tensors) == num_expect
