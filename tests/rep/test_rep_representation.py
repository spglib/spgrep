from __future__ import annotations

import numpy as np

from spgrep.rep.irreps import is_equivalent_irrep
from spgrep.rep.representation import get_character, get_intertwiner


def test_intertwiner():
    rep1 = np.array(
        [
            [[1.0 - 0.0j, 0.0 + 0.0j], [0.0 + 0.0j, 1.0 - 0.0j]],
            [[0.0 + 0.0j, 1.0 - 0.0j], [1.0 - 0.0j, 0.0 + 0.0j]],
            [[-0.0 - 1.0j, 0.0 + 0.0j], [0.0 + 0.0j, -0.0 + 1.0j]],
            [[0.0 + 0.0j, -0.0 - 1.0j], [-0.0 + 1.0j, 0.0 + 0.0j]],
        ]
    )
    rep2 = np.array(
        [
            [[1.0 - 0.0j, 0.0 + 0.0j], [0.0 + 0.0j, 1.0 - 0.0j]],
            [[0.0 + 0.0j, 1.0 - 0.0j], [1.0 - 0.0j, 0.0 + 0.0j]],
            [[-0.0 + 1.0j, 0.0 - 0.0j], [0.0 - 0.0j, 0.0 - 1.0j]],
            [[0.0 - 0.0j, -0.0 + 1.0j], [0.0 - 1.0j, 0.0 - 0.0j]],
        ]
    )

    intertwiner = get_intertwiner(rep1, rep2)
    assert is_equivalent_irrep(get_character(rep1), get_character(rep2))
    assert np.allclose(
        np.einsum("kil,lj->kij", rep1, intertwiner),
        np.einsum("il,klj->kij", intertwiner, rep2),
    )
