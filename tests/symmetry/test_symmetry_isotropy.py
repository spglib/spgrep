from __future__ import annotations

import numpy as np

from spgrep.symmetry.group import get_cayley_table
from spgrep.symmetry.pointgroup import pg_dataset
from spgrep.symmetry.subgroup import enumerate_point_subgroup, enumerate_point_subgroup_naive


def test_enumerate_point_subgroup():
    pointgroup = pg_dataset["4/mmm"][0]
    table = get_cayley_table(np.array(pointgroup))
    flags = [True for _ in range(len(pointgroup))]
    subgroups_actual = enumerate_point_subgroup(table, flags, return_conjugacy_class=False)
    subgroups_expect = enumerate_point_subgroup_naive(table, flags)
    assert len(subgroups_actual) == len(subgroups_expect)
