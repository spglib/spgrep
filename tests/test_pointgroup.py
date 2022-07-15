import numpy as np

from spgrep.group import get_cayley_table, is_matrix_group
from spgrep.pointgroup import pg_dataset


def test_pg_dataset():
    for _, groups in pg_dataset.items():
        for pg in groups:
            rotations = np.array(pg)
            assert is_matrix_group(rotations)
            table = get_cayley_table(rotations)
            assert table.shape[0] == len(rotations)
