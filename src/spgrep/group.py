"""Group-theory related functions."""

from __future__ import annotations

import warnings

warnings.warn(
    "This module is deprecated and will be removed in a future release. "
    "Use `spgrep.rep.group` (finite-group utilities) and/or "
    "`spgrep.symmetry.group` (space-group/little-group utilities) instead.",
    DeprecationWarning,
    stacklevel=2,
)

from spgrep.rep.group import get_identity_index, get_inverse_index, get_order  # noqa: E402
from spgrep.symmetry.group import (  # noqa: E402
    check_cocycle_condition,
    decompose_by_maximal_space_subgroup,
    get_cayley_table,
    get_factor_system_from_little_group,
    get_little_group,
    is_matrix_group,
)

__all__ = [
    "get_identity_index",
    "get_inverse_index",
    "get_order",
    "get_cayley_table",
    "is_matrix_group",
    "get_factor_system_from_little_group",
    "get_little_group",
    "check_cocycle_condition",
    "decompose_by_maximal_space_subgroup",
]
