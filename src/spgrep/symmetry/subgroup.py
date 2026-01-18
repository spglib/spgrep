"""Isotropy subgroup of space group."""

from __future__ import annotations

from queue import Queue

from spgrep.rep.group import get_identity_index, get_inverse_index
from spgrep.utils import (
    NDArrayInt,
)


def enumerate_point_subgroup(
    table: NDArrayInt, preserve_sublattice: list[bool], return_conjugacy_class: bool = True
) -> list[list[int]]:
    """Enumerate conjugacy subgroups of point group.

    Parameters
    ----------
    table: array[int], (order, order)
        Multiplication table of group
    preserve_sublattice: list[bool]
        Specify ``preserve_sublattice[i] = True`` if the ``i``-th operation preserves translational subgroup of isotropy subgroup
    return_conjugacy_class: bool, default=True
        If true, return representatives of conjugacy classes.

    Returns
    -------
    subgroups: list[list[int]]
    """
    order = len(table)
    identity = get_identity_index(table)
    # Represent choice of elements by bit array
    st = {1 << identity}
    for i in range(order):
        if (not preserve_sublattice[i]) or (i == identity):
            continue
        if (1 << i) in st:
            # Already visited
            continue

        next_st = set()
        for bits in st:
            elements = _decode_bits(bits, order)
            assert _is_subgroup(elements, table)
            generated = _traverse(elements + [i], identity, table)
            next_st.add(sum(1 << idx for idx in set(generated)))

        st = st.union(next_st)

    if not return_conjugacy_class:
        subgroups = []
        for bits in sorted(st):
            subgroups.append(_decode_bits(bits, order))
        return subgroups

    # Group by conjugacy classes
    found = set()
    ret = []
    for bits in sorted(st):
        if bits in found:
            continue
        elements = _decode_bits(bits, order)
        ret.append(elements)
        for i in range(order):
            if not preserve_sublattice[i]:
                continue
            inv = get_inverse_index(table, i)
            conj = [int(table[inv, table[idx, i]]) for idx in elements]
            found.add(sum(1 << idx for idx in set(conj)))

    assert found == st
    return ret


def enumerate_point_subgroup_naive(table, preserve_sublattice: list[bool]):
    """Enumerate conjugacy subgroups of point group in brute force."""
    order = len(table)
    ret = []
    for bits in range(1, 1 << order):
        elements = _decode_bits(bits, order)
        if not all([preserve_sublattice[idx] for idx in elements]):
            continue

        if _is_subgroup(elements, table):
            ret.append(bits)

    return ret


def _decode_bits(bits: int, order: int) -> list[int]:
    elements = [idx for idx in range(order) if (bits >> idx) & 1 == 1]
    return elements


def _is_subgroup(elements: list[int], table: NDArrayInt) -> bool:
    subtable = table[elements][:, elements]
    for i in range(len(subtable)):
        if (set(subtable[i]) != set(elements)) or (set(subtable[:, i]) != set(elements)):
            return False
    return True


def _traverse(
    generators: list[int],
    identity: int,
    table: NDArrayInt,
) -> list[int]:
    """Traverse group elements from generators."""
    visited = [False for _ in range(len(table))]
    que = Queue()  # type: ignore
    que.put(identity)

    while not que.empty():
        g = que.get()
        if visited[g]:
            continue
        visited[g] = True

        for h in generators:
            gh = int(table[g, h])  # cast np.int64 to int
            if not visited[gh]:
                que.put(gh)

    return sorted([i for i, v in enumerate(visited) if v])
