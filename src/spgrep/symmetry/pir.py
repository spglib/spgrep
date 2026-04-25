"""Physically irreducible representations of space-group little groups."""

from __future__ import annotations

import numpy as np

from spgrep._constants import ATOL, MAX_NUM_RANDOM_GENERATIONS
from spgrep.rep.irreps import frobenius_schur_indicator, is_equivalent_irrep
from spgrep.rep.pir import get_physically_irrep
from spgrep.rep.representation import get_character
from spgrep.utils import NDArrayBool, NDArrayComplex, NDArrayFloat


def purify_real_irrep_value(real_irrep: NDArrayFloat, atol: float = ATOL) -> NDArrayFloat:
    """Purify values of physically irreducible representations."""
    values = [
        0,
        1,  # 0/1
        1 / 2,  # 1/3
        np.sqrt(3) / 2,  # 1/3
        -1 / 2,  # 2/3
        -np.sqrt(3) / 2,  # 2/3
        -1,  # 1/2
    ]
    for v in values:
        real_irrep[np.abs(real_irrep - v) < atol] = v
    return real_irrep


def _make_physically_irreps_from_complex(
    irreps: list[NDArrayComplex],
    kpoint: NDArrayFloat | None = None,
    atol: float = ATOL,
    max_num_random_generations: int = MAX_NUM_RANDOM_GENERATIONS,
) -> tuple[list[NDArrayFloat], list[int]]:
    """Convert a list of complex unitary irreps into physically irreducible real forms.

    Pairs conjugate irreps, computes the Frobenius-Schur indicator of each, and
    realifies via :func:`spgrep.rep.pir.get_physically_irrep`. When ``kpoint`` is
    given the caller is operating on a small representation of a little group
    and the character-sum formula for the indicator vanishes unless
    ``2 * kpoint`` is an integer vector; in that case the indicator is forced to
    0 (Re/Im block doubling).
    """
    conjugated_pairs = []
    visited = [False for _ in range(len(irreps))]
    characters = [get_character(irrep) for irrep in irreps]
    for i, ci in enumerate(characters):
        if visited[i]:
            continue
        visited[i] = True
        inequivalent = False
        for j, cj in enumerate(characters):
            if visited[j]:
                continue
            if is_equivalent_irrep(np.conj(ci), cj):
                conjugated_pairs.append((i, j))
                visited[j] = True
                inequivalent = True
                break
        if not inequivalent:
            conjugated_pairs.append((i, i))

    two_kpoint_is_zero = True
    if kpoint is not None:
        two_kpoint = 2 * kpoint
        two_kpoint -= np.rint(two_kpoint)
        two_kpoint_is_zero = bool(np.allclose(two_kpoint, 0))

    real_irreps: list[NDArrayFloat] = []
    indicators: list[int] = []
    for conj_pair in conjugated_pairs:
        irrep = irreps[conj_pair[0]]
        indicator = frobenius_schur_indicator(irrep) if two_kpoint_is_zero else 0

        real_irrep = get_physically_irrep(
            irrep, indicator, atol=atol, max_num_random_generations=max_num_random_generations
        )
        real_irrep = purify_real_irrep_value(real_irrep, atol=atol)
        real_irreps.append(real_irrep)
        indicators.append(indicator)

    return real_irreps, indicators


def _realify_corep_to_real_rep(
    corep: NDArrayComplex,
    anti_linear: NDArrayBool,
) -> NDArrayFloat:
    r"""Realise a complex corepresentation as a real linear matrix representation.

    A corepresentation acts on a complex vector space :math:`\mathbb{C}^d` via
    unitary matrices on linear ops and via :math:`z \mapsto D(g) z^{\ast}` on
    antilinear ops. Viewing :math:`\mathbb{C}^d \simeq \mathbb{R}^{2d}` via
    :math:`z = x + i y \mapsto (x, y)`, multiplication by :math:`w = a + i b` is
    :math:`\bigl[\begin{smallmatrix}a & -b \\ b & a\end{smallmatrix}\bigr]` and
    complex conjugation is :math:`\mathrm{diag}(I_d, -I_d)`. Composing these
    gives the real-linear images:

    - linear op :math:`g`:
      :math:`\bigl[\begin{smallmatrix} \mathrm{Re}\,D & -\mathrm{Im}\,D \\
      \mathrm{Im}\,D & \mathrm{Re}\,D \end{smallmatrix}\bigr]`,
    - antilinear op :math:`g`:
      :math:`\bigl[\begin{smallmatrix} \mathrm{Re}\,D & \mathrm{Im}\,D \\
      \mathrm{Im}\,D & -\mathrm{Re}\,D \end{smallmatrix}\bigr]`.

    This conversion produces a genuine matrix representation of the underlying
    (linear) group on :math:`\mathbb{R}^{2d}`, which is the physically
    irreducible representation associated with the corepresentation.

    Parameters
    ----------
    corep: array, (order, dim, dim)
        Complex corepresentation matrices.
    anti_linear: array[bool], (order,)
        ``anti_linear[g]`` is ``True`` iff operation ``g`` acts antilinearly
        (i.e. its corepresentation product law involves complex conjugation).

    Returns
    -------
    real_rep: array, (order, 2 * dim, 2 * dim)
    """
    order, dim, _ = corep.shape
    real_rep = np.empty((order, 2 * dim, 2 * dim), dtype=np.float64)
    re = np.real(corep)
    im = np.imag(corep)
    for g in range(order):
        if anti_linear[g]:
            real_rep[g, :dim, :dim] = re[g]
            real_rep[g, :dim, dim:] = im[g]
            real_rep[g, dim:, :dim] = im[g]
            real_rep[g, dim:, dim:] = -re[g]
        else:
            real_rep[g, :dim, :dim] = re[g]
            real_rep[g, :dim, dim:] = -im[g]
            real_rep[g, dim:, :dim] = im[g]
            real_rep[g, dim:, dim:] = re[g]
    return real_rep
