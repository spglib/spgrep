from __future__ import annotations

import warnings

warnings.warn(
    "This module is deprecated and will be removed in a future release. "
    "Use `spgrep.symmetry.pointgroup` instead.",
    DeprecationWarning,
    stacklevel=2,
)

from spgrep.symmetry.pointgroup import (  # noqa: E402
    get_generators,
    get_pointgroup_chain_generators,
    pg_dataset,
    pg_solvable_chain,
)

__all__ = [
    "pg_dataset",
    "pg_solvable_chain",
    "get_generators",
    "get_pointgroup_chain_generators",
]
