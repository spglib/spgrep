from __future__ import annotations

import numpy as np
import pytest

from spgrep.pointgroup import pg_dataset
from spgrep.utils import NDArrayFloat, NDArrayInt, get_symmetry_from_hall_number


@pytest.fixture
def Oh() -> NDArrayInt:
    rotations = np.array(pg_dataset["m-3m"][0])
    return rotations


@pytest.fixture
def C4() -> NDArrayInt:
    rotations = np.array(pg_dataset["4"][0])
    return rotations


@pytest.fixture
def C3v() -> NDArrayInt:
    # g0
    # g1 = 3+
    # g2 = g1^-1
    # g3 = m
    # g4 = g1^-1 * g3
    # g5 = g1 * g3
    rotations = np.array(pg_dataset["3m"][0])  # 3m1
    return rotations


@pytest.fixture
def P42mnm() -> tuple[NDArrayInt, NDArrayFloat]:
    # P4_2/mnm (No. 136)
    return get_symmetry_from_hall_number(hall_number=419)


@pytest.fixture
def Ia3d() -> tuple[NDArrayInt, NDArrayFloat]:
    # Ia3d (No. 230)
    return get_symmetry_from_hall_number(hall_number=530)


@pytest.fixture
def corundum_cell():
    # Corundum structure, R-3c (No. 167)
    # https://materialsproject.org/materials/mp-1143
    a = 4.80502783
    c = 13.11625361
    lattice = np.array(
        [
            [a, 0, 0],
            [-1 / 2 * a, np.sqrt(3) / 2 * a, 0],
            [0, 0, c],
        ]
    )
    positions = np.array(
        [
            [0.00000000, 0.00000000, 0.14790400],
            [0.33333333, 0.66666667, 0.01876267],
            [0.33333333, 0.66666667, 0.31457067],
            [0.66666667, 0.33333333, 0.18542933],
            [0.66666667, 0.33333333, 0.48123733],
            [0.00000000, 0.00000000, 0.35209600],
            [0.00000000, 0.00000000, 0.64790400],
            [0.33333333, 0.66666667, 0.51876267],
            [0.33333333, 0.66666667, 0.81457067],
            [0.66666667, 0.33333333, 0.68542933],
            [0.66666667, 0.33333333, 0.98123733],
            [0.00000000, 0.00000000, 0.85209600],
            [0.30614600, 0.00000000, 0.25000000],
            [0.66666667, 0.02718733, 0.08333333],
            [0.00000000, 0.30614600, 0.25000000],
            [0.69385400, 0.69385400, 0.25000000],
            [0.97281267, 0.63947933, 0.08333333],
            [0.36052067, 0.33333333, 0.08333333],
            [0.97281267, 0.33333333, 0.58333333],
            [0.33333333, 0.36052067, 0.41666667],
            [0.66666667, 0.63947933, 0.58333333],
            [0.36052067, 0.02718733, 0.58333333],
            [0.63947933, 0.97281267, 0.41666667],
            [0.02718733, 0.66666667, 0.41666667],
            [0.63947933, 0.66666667, 0.91666667],
            [0.00000000, 0.69385400, 0.75000000],
            [0.33333333, 0.97281267, 0.91666667],
            [0.02718733, 0.36052067, 0.91666667],
            [0.30614600, 0.30614600, 0.75000000],
            [0.69385400, 0.00000000, 0.75000000],
        ]
    )
    numbers = [0] * 12 + [1] * 18

    return (lattice, positions, numbers)
