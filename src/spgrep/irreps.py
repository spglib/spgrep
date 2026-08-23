"""On-the-fly irreps generations."""

from __future__ import annotations

import warnings

warnings.warn(
    "spgrep.irreps is deprecated and will be removed in future releases. "
    "Please use spgrep.rep.enumerate, spgrep.rep.pir, and spgrep.symmetry.enumerate instead.",
    DeprecationWarning,
    stacklevel=2,
)

from spgrep.rep.enumerate import (  # noqa: E402
    decompose_representation,
    enumerate_unitary_irreps_from_regular_representation,
)
from spgrep.rep.pir import get_physically_irrep  # noqa: E402
from spgrep.symmetry.enumerate import (  # noqa: E402
    enumerate_small_representations,
    enumerate_unitary_irreps,
    enumerate_unitary_irreps_from_solvable_group_chain,
    purify_irrep_value,
)
from spgrep.symmetry.pir import purify_real_irrep_value  # noqa: E402

__all__ = [
    "decompose_representation",
    "enumerate_small_representations",
    "enumerate_unitary_irreps",
    "enumerate_unitary_irreps_from_regular_representation",
    "enumerate_unitary_irreps_from_solvable_group_chain",
    "get_physically_irrep",
    "purify_irrep_value",
    "purify_real_irrep_value",
]
