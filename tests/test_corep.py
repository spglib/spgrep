import numpy as np
import pytest

from spgrep.core import (
    get_crystallographic_pointgroup_spinor_irreps_from_symmetry,
    get_magnetic_spacegroup_coreps_from_primitive_symmetry,
)
from spgrep.corep import (
    _enumerate_small_corepresentations_with_factor_system,
    get_corep_spinor_factor_system,
)
from spgrep.symmetry.group import get_cayley_table, get_factor_system_from_little_group
from spgrep.utils import (
    NDArrayBool,
    NDArrayComplex,
    NDArrayFloat,
    NDArrayInt,
    get_symmetry_from_hall_number,
)


def _assert_corep_composition_law(
    corep: NDArrayComplex,
    anti_linear: NDArrayBool,
    factor_system: NDArrayComplex,
    table: NDArrayInt,
    atol: float = 1e-8,
) -> None:
    """Verify ``D(g_i) @ (D(g_j)* if g_i anti-linear else D(g_j)) == omega(g_i, g_j) * D(g_i g_j)``."""
    order = corep.shape[0]
    for i in range(order):
        for j in range(order):
            right = np.conj(corep[j]) if anti_linear[i] else corep[j]
            lhs = corep[i] @ right
            rhs = factor_system[i, j] * corep[table[i, j]]
            assert np.allclose(lhs, rhs, atol=atol), f"corep law violated at (i={i}, j={j})"


def check_corep_cocycle_condition(
    rotations: NDArrayInt,
    time_reversals: NDArrayInt,
    corep_factor_system: NDArrayComplex,
) -> bool:
    """Return true if given factor system satisfies the cocycle condition for co-representation."""
    order = len(rotations)
    table = get_cayley_table(rotations, time_reversals)

    # Check associativity: (gi * gj) * gk = gi * (gj * gk)
    for i, tri in enumerate(time_reversals):
        for j in range(order):
            for k in range(order):
                ij = table[i, j]
                jk = table[j, k]

                if tri == 1:
                    if not np.isclose(
                        corep_factor_system[i, j] * corep_factor_system[ij, k],
                        corep_factor_system[i, jk] * np.conj(corep_factor_system[j, k]),
                    ):
                        return False
                else:
                    if not np.isclose(
                        corep_factor_system[i, j] * corep_factor_system[ij, k],
                        corep_factor_system[i, jk] * corep_factor_system[j, k],
                    ):
                        return False

    return True


@pytest.mark.parametrize(
    "symmetry_and_lattice",
    [
        ("P42mnm_type3"),
        # ("bcc_type4"),  # Passed but too long to test
    ],
)
def test_corep_spinor_factor_system(request, symmetry_and_lattice):
    rotations, _, time_reversals, lattice = request.getfixturevalue(symmetry_and_lattice)

    corep_spinor_factor_system, unitary_rotations, anti_linear = get_corep_spinor_factor_system(
        lattice, rotations, time_reversals
    )

    # Check unitary
    for unitary_rotation in unitary_rotations:
        assert np.allclose(
            unitary_rotation @ np.conj(unitary_rotation).T,
            np.eye(2, dtype=np.complex128),
        )

    # Cocycle condition
    assert check_corep_cocycle_condition(rotations, time_reversals, corep_spinor_factor_system)


@pytest.mark.parametrize("method", [("Neto"), ("random")])
@pytest.mark.parametrize(
    "symmetry_and_lattice",
    [
        ("P42mnm_type1"),
        ("P42mnm_type2"),
        ("P42mnm_type3"),
        ("bcc_type4"),
    ],
)
def test_get_crystallographic_pointgroup_spinor_irreps_from_symmetry(
    request, method, symmetry_and_lattice
):
    rotations, _, time_reversals, lattice = request.getfixturevalue(symmetry_and_lattice)

    # TODO: Add more tests
    (
        co_irreps,
        indicators,
        factor_system,
        unitary_rotations,
        anti_linear,
    ) = get_crystallographic_pointgroup_spinor_irreps_from_symmetry(
        lattice,
        rotations,
        time_reversals,
        method=method,
    )

    table = get_cayley_table(rotations, time_reversals)
    for co_irrep in co_irreps:
        _assert_corep_composition_law(co_irrep, anti_linear, factor_system, table)


@pytest.mark.parametrize(
    "kpoint",
    [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.5, 0.0, 0.0]),
        np.array([0.5, 0.5, 0.0]),
        np.array([0.5, 0.5, 0.5]),
    ],
)
@pytest.mark.parametrize("method", [("Neto"), ("random")])
@pytest.mark.parametrize(
    "symmetry_and_lattice",
    [
        ("P42mnm_type1"),
        ("P42mnm_type2"),
        ("P42mnm_type3"),
        ("bcc_type4"),
    ],
)
def test_get_magnetic_spacegroup_coreps_from_primitive_symmetry(
    request, kpoint: NDArrayFloat, method, symmetry_and_lattice
):
    rotations, translations, time_reversals, lattice = request.getfixturevalue(
        symmetry_and_lattice
    )

    # TODO: Add more tests
    (
        co_irreps,
        anti_linear,
        mapping_little_group,
    ) = get_magnetic_spacegroup_coreps_from_primitive_symmetry(
        rotations=rotations,
        translations=translations,
        time_reversals=time_reversals,
        kpoint=kpoint,
        method=method,
    )


def test_corep_composition_law_pm3m_lambda():
    """Regression: co-rep extension must conjugate the right operand when the left is anti-linear.

    Constructs the little group of :math:`\\pm \\mathbf{k}` for Pm-3m at
    :math:`\\mathbf{k} = (1/4, 1/4, 1/4)` (a non-TRIM point; the 2-dim ``E`` irrep
    of :math:`D_3` appears here and is non-real). Feeds the k-flip mask as
    ``little_time_reversals`` -- this matches how the co-rep machinery is used to
    extend a small rep of :math:`G^{\\mathbf{k}}` to :math:`G^{\\mathbf{k}\\bar{\\mathbf{k}}}`.
    Before the conjugation fix this assertion fails at 72/144 pairs for one of
    the three ``ksi`` blocks; after the fix all pairs pass.
    """
    rotations, translations = get_symmetry_from_hall_number(517)  # Pm-3m
    kpoint = np.array([0.25, 0.25, 0.25])

    pm_k_mask = []
    flip_k_list = []
    for rot in rotations:
        residual_plus = rot.T @ kpoint - kpoint
        residual_minus = rot.T @ kpoint + kpoint
        keeps_k = np.allclose(residual_plus - np.rint(residual_plus), 0)
        flips_k = np.allclose(residual_minus - np.rint(residual_minus), 0)
        if keeps_k:
            pm_k_mask.append(True)
            flip_k_list.append(0)
        elif flips_k:
            pm_k_mask.append(True)
            flip_k_list.append(1)
        else:
            pm_k_mask.append(False)
    pm_k_mask = np.array(pm_k_mask)
    flip_k = np.array(flip_k_list, dtype=int)
    pm_k_rotations = rotations[pm_k_mask]
    pm_k_translations = translations[pm_k_mask]

    factor_system = get_factor_system_from_little_group(pm_k_rotations, pm_k_translations, kpoint)
    table = get_cayley_table(pm_k_rotations, flip_k)

    coreps, _ = _enumerate_small_corepresentations_with_factor_system(
        little_rotations=pm_k_rotations,
        little_translations=pm_k_translations,
        little_time_reversals=flip_k,
        kpoint=kpoint,
        factor_system=factor_system,
    )

    # Translation phases carried on the small co-reps -- strip them so the
    # residual composition law uses the linear-part factor system directly.
    phases = np.exp(-2j * np.pi * pm_k_translations @ kpoint)
    anti_linear = flip_k.astype(bool)
    for corep in coreps:
        bare = corep / phases[:, None, None]
        _assert_corep_composition_law(bare, anti_linear, factor_system, table)
