from __future__ import annotations

import warnings

warnings.warn(
    "The module 'spgrep.transform' is deprecated and will be removed in future versions. "
    "Please use 'spgrep.symmetry.transform' instead.",
    DeprecationWarning,
    stacklevel=2,
)

from spgrep.symmetry.transform import (  # noqa: E402
    get_centering,
    get_crystal_system,
    get_primitive_transformation_matrix,
    transform_symmetry_and_kpoint,
    unique_primitive_symmetry,
)

__all__ = [
    "transform_symmetry_and_kpoint",
    "unique_primitive_symmetry",
    "get_primitive_transformation_matrix",
    "get_crystal_system",
    "get_centering",
]
