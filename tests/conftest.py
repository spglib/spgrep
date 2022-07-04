import numpy as np
import pytest

from spgrep.utils import NDArrayInt


@pytest.fixture
def C3v() -> NDArrayInt:
    # fmt: off
    rotations = np.array([
        [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ],
        [
            [0, -1, 0],
            [1, -1, 0],
            [0, 0, 1],
        ],
        [
            [-1, 1, 0],
            [-1, 0, 0],
            [0, 0, 1],
        ],
        [
            [0, -1, 0],
            [-1, 0, 0],
            [0, 0, 1],
        ],
        [
            [-1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
        ],
        [
            [1, 0, 0],
            [1, -1, 0],
            [0, 0, 1],
        ],
    ])
    # fmt: on
    return rotations
