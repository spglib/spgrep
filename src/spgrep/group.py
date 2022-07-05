from __future__ import annotations

from itertools import product

import numpy as np

from spgrep.utils import NDArrayInt, is_integer_array, ndarray2d_to_integer_tuple


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
