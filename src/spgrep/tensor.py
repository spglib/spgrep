from __future__ import annotations

import warnings

from spgrep.symmetry.tensor import (
    apply_intrinsic_symmetry,
    get_symmetry_adapted_tensors,
)

warnings.warn(
    "spgrep.symmetry.tensor is deprecated. Please use spgrep.symmetry.tensor instead.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = [
    "apply_intrinsic_symmetry",
    "get_symmetry_adapted_tensors",
]
