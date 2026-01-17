"""Group-theory related functions."""

from __future__ import annotations

import numpy as np

from spgrep.utils import (
    NDArrayInt,
)


def get_identity_index(table: NDArrayInt) -> int:
    """Return index for identity of group."""
    order = table.shape[0]
    for i in range(order):
        if np.all(table[i, :] == np.arange(order)):
            return i

    raise ValueError("Unreachable!")


def get_inverse_index(table: NDArrayInt, idx: int) -> int:
    """Return index of inverse of ``idx`` element in ``table``."""
    order = table.shape[0]
    id_idx = get_identity_index(table)
    for i in range(order):
        if table[idx, i] == id_idx:
            return i

    raise ValueError("Unreachable!")


def get_order(table: NDArrayInt, idx: int) -> int:
    """Return order of element ``idx`` in ``table``."""
    id_idx = get_identity_index(table)
    ret = 1
    tmp = idx
    while tmp != id_idx:
        tmp = table[tmp, idx]
        ret += 1
    return ret
