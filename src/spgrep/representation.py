"""Representation-matrix related implementations."""

from __future__ import annotations

import warnings

warnings.warn(
    "The `spgrep.representation` module is deprecated and will be removed in a future release. "
    "Use `spgrep.rep.representation`, `spgrep.symmetry.representation`, and `spgrep.rep.irreps` instead.",
    DeprecationWarning,
    stacklevel=2,
)


from spgrep.rep.irreps import frobenius_schur_indicator, project_to_irrep  # noqa: E402
from spgrep.rep.representation import (  # noqa: E402
    get_character,
    get_direct_product,
    get_intertwiner,
    is_representation,
    is_unitary,
)
from spgrep.symmetry.representation import (  # noqa: E402
    check_spacegroup_representation,
    get_projective_regular_representation,
    get_regular_representation,
)

__all__ = [
    # spgrep.rep.representation
    "get_intertwiner",
    "get_character",
    "is_unitary",
    "is_representation",
    "get_direct_product",
    # spgrep.symmetry.representation
    "get_regular_representation",
    "get_projective_regular_representation",
    "check_spacegroup_representation",
    # spgrep.rep.irreps
    "project_to_irrep",
    "frobenius_schur_indicator",
]
