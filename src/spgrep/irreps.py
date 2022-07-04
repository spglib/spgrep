import numpy as np

from spgrep.utils import NDArrayInt, ndarray2d_to_integer_tuple


def get_regular_representation(rotations: NDArrayInt) -> NDArrayInt:
    """Calculate regular representation of point group.

    Parameters
    ----------
    rotations: array, (num_sym, 3, 3)

    Returns
    -------
    reg: array, (num_sym, num_sym, num_sym)
        ``reg[k]`` is a representation matrix for ``rotations[k]``.
        If and only if ``np.dot(rotations[k], rotations[j]) == rotations[i]``, ``reg[k, i, j] == 1``.
    """
    n = len(rotations)
    reg = np.zeros((n, n, n), dtype=int)
    rotations_list = [ndarray2d_to_integer_tuple(r) for r in rotations]

    for k, gk in enumerate(rotations):
        for j, gj in enumerate(rotations):
            gkj = np.dot(gk, gj)
            try:
                i = rotations_list.index(ndarray2d_to_integer_tuple(gkj))
            except ValueError:
                raise ValueError("Given matrices should form group.")

            reg[k, j] = i

    return reg
