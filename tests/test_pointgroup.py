import numpy as np

from spgrep.group import get_cayley_table, get_inverse_index, is_matrix_group
from spgrep.pointgroup import get_generators, pg_dataset


def test_pg_dataset():
    for _, groups in pg_dataset.items():
        for pg in groups:
            rotations = np.array(pg)
            assert is_matrix_group(rotations)
            table = get_cayley_table(rotations)
            assert table.shape[0] == len(rotations)


def test_group_chain():
    for pg_symbol, groups in pg_dataset.items():
        for idx, pg in enumerate(groups):
            gens = get_generators(pg_symbol, idx)

            order = len(pg)
            table = get_cayley_table(np.array(pg))
            inv = [get_inverse_index(table, g) for g in range(order)]

            prev_size = order
            for gen in gens:
                # table[:gen, :gen] is normal subgroup
                assert np.all(table[:gen, :gen] < gen)
                for n in range(gen):
                    assert all([table[inv[g], table[n, g]] < gen for g in range(prev_size)])
                prev_size = gen
