import numpy as np
import pytest

from spgrep.utils import NDArrayInt


@pytest.fixture
def C3v() -> NDArrayInt:
    # fmt: off
    rotations = np.array([
        # g0
        [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ],
        # g1 = 3+
        [
            [0, -1, 0],
            [1, -1, 0],
            [0, 0, 1],
        ],
        # g2 = g1^-1
        [
            [-1, 1, 0],
            [-1, 0, 0],
            [0, 0, 1],
        ],
        # g3 = m
        [
            [0, -1, 0],
            [-1, 0, 0],
            [0, 0, 1],
        ],
        # g4 = g1^-1 * g3
        [
            [-1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
        ],
        # g5 = g1 * g3
        [
            [1, 0, 0],
            [1, -1, 0],
            [0, 0, 1],
        ],
    ])
    # fmt: on
    return rotations
