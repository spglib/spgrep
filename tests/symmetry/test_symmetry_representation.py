from __future__ import annotations

import numpy as np

from spgrep.rep.representation import get_character
from spgrep.symmetry.representation import get_regular_representation


def test_get_character(C3v):
    reg = get_regular_representation(C3v)
    actual = get_character(reg)
    expect = np.array([6, 0, 0, 0, 0, 0])
    assert np.allclose(actual, expect)
