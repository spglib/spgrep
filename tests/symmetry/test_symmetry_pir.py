import numpy as np
import pytest

from spgrep.core import get_spacegroup_irreps_from_primitive_symmetry
from spgrep.rep.representation import is_representation
from spgrep.symmetry.enumerate import enumerate_small_representations, enumerate_unitary_irreps
from spgrep.symmetry.group import (
    get_cayley_table,
    get_little_group,
    get_little_group_of_pm_k,
)
from spgrep.utils import (
    NDArrayComplex,
    NDArrayFloat,
    NDArrayInt,
    get_symmetry_from_hall_number,
    ndarray2d_to_integer_tuple,
)


@pytest.mark.parametrize(
    "fixture_name,num_expect",
    [
        ("C4", 3),
        ("Oh", 10),
    ],
)
def test_physically_irreducible_representation(request, fixture_name, num_expect):
    rotations = request.getfixturevalue(fixture_name)
    real_irreps, _ = enumerate_unitary_irreps(rotations, real=True)
    assert len(real_irreps) == num_expect

    # Check representation's property
    table = get_cayley_table(rotations)
    for irrep in real_irreps:
        assert is_representation(irrep, table)


def check_spacegroup_pir(
    little_rotations: NDArrayInt,
    little_translations: NDArrayFloat,
    kpoint: NDArrayFloat,
    rep: NDArrayComplex,
    indicator: int,
    rtol: float = 1e-5,
) -> bool:
    little_rotations_int = [ndarray2d_to_integer_tuple(rotation) for rotation in little_rotations]

    # Check if ``rep`` preserves multiplication
    for r1, t1, m1 in zip(little_rotations, little_translations, rep):
        for r2, t2, m2 in zip(little_rotations, little_translations, rep):
            # (r1, t1) (r2, t2) = (E, additional) (r12, little_translations[idx])
            r12 = r1 @ r2
            t12 = r1 @ t2 + t1

            # little_translations[idx] may differ from t12 by lattice translation.
            idx = little_rotations_int.index(ndarray2d_to_integer_tuple(r12))
            additional = t12 - little_translations[idx]
            if indicator == 1:
                # Must be 1 or -1
                phase = np.real(np.exp(-2j * np.pi * np.dot(kpoint, additional)))
                m12 = phase * rep[idx]
            else:
                half_dim = rep.shape[1] // 2
                rep_additional = np.zeros((2 * half_dim, 2 * half_dim), dtype=np.complex128)
                # [ [ cos(kt), -sin(kt) ],
                #   [ sin(kt), cos(kt) ] ]
                rep_additional[:half_dim, :half_dim] = np.cos(
                    2 * np.pi * np.dot(kpoint, additional)
                ) * np.eye(half_dim)
                rep_additional[half_dim:, :half_dim] = np.sin(
                    2 * np.pi * np.dot(kpoint, additional)
                ) * np.eye(half_dim)
                rep_additional[half_dim:, half_dim:] = rep_additional[:half_dim, :half_dim]
                rep_additional[:half_dim:, half_dim:] = -rep_additional[half_dim:, :half_dim]

                m12 = rep_additional @ rep[idx]

            if not np.allclose(m12, m1 @ m2, rtol=rtol):
                return False

    return True


def test_spacegroup_pir():
    # P6_1 (No. 169)
    rotations, translations = get_symmetry_from_hall_number(463)
    kpoint = np.array([1 / 3, 1 / 3, 1 / 2])  # H point
    little_rotations, little_translations, _ = get_little_group(rotations, translations, kpoint)
    assert len(little_rotations) == 3

    physically_irreps, indicators = enumerate_small_representations(
        little_rotations, little_translations, kpoint, real=True
    )
    assert len(physically_irreps) == 2
    for pir, indicator in zip(physically_irreps, indicators):
        assert check_spacegroup_pir(little_rotations, little_translations, kpoint, pir, indicator)


def test_spacegroup_pir_pm3m_lambda():
    # SG 221 (Pm-3m) at k=(u,u,u) on the Lambda line, where 2k is not
    # equivalent to 0. Compares against the Bilbao Crystallographic Server's
    # PIR enumeration (LD1, LD2, LD3 of dimensions 2, 2, 4 on the little group
    # of +/- k). We verify the structural properties that must hold for any
    # valid PIR extension on G^{k,k-bar}; the exact matrix entries are
    # basis-dependent.
    rotations, translations = get_symmetry_from_hall_number(517)  # Pm-3m
    u = 0.1  # generic Lambda point; 2u is not integer so 2k is not equivalent to 0
    kpoint = np.array([u, u, u])

    pm_k_rotations, _, mapping_little_group, flip_k = get_little_group_of_pm_k(
        rotations, translations, kpoint
    )
    # little cogroup of +/- k for Lambda in Pm-3m is D3d (order 12),
    # half stabilizing k, half mapping k to -k.
    assert len(pm_k_rotations) == 12
    assert int(np.sum(flip_k == 0)) == 6
    assert int(np.sum(flip_k == 1)) == 6

    physically_irreps, _ = get_spacegroup_irreps_from_primitive_symmetry(
        rotations, translations, kpoint, real=True
    )
    # Bilbao reports LD1, LD2, LD3 -- 3 PIRs in total.
    assert len(physically_irreps) == 3

    table = get_cayley_table(pm_k_rotations)
    for pir in physically_irreps:
        # Each output must be a real matrix representation of G^{k,k-bar}.
        assert pir.dtype.kind == "f"
        assert pir.shape[0] == len(pm_k_rotations)
        assert pir.shape[1] == pir.shape[2]
        assert is_representation(pir, table)


def test_spacegroup_pir_r3m_smk1():
    # SG 166 (R-3m) at Bilbao's k1=(u,-2u,0)_hex on the SM line, where 2k is
    # not equivalent to 0. Compares against the Bilbao Crystallographic Server's
    # PIR enumeration (SM1, SM2 of dimension 2 on the 4-op little group of
    # +/- k). As in test_spacegroup_pir_pm3m_lambda, only the structural
    # properties are checked; the exact matrix entries are basis-dependent.
    # Hall 459 is the primitive rhombohedral setting; Bilbao's hexagonal-axes
    # k1=(u,-2u,0) transforms to k=(0,-u,u) in primitive rhombohedral basis.
    rotations, translations = get_symmetry_from_hall_number(459)  # R-3m primitive rhom
    u = 0.13  # generic; 2u is not an integer so 2k is not equivalent to 0
    kpoint = np.array([0.0, -u, u])

    pm_k_rotations, _, _, flip_k = get_little_group_of_pm_k(rotations, translations, kpoint)
    # Little group of +/- k1 = {1, 2_010, 1bar, m_010}: 2 stabilize k1, 2 flip.
    assert len(pm_k_rotations) == 4
    assert int(np.sum(flip_k == 0)) == 2
    assert int(np.sum(flip_k == 1)) == 2

    physically_irreps, _ = get_spacegroup_irreps_from_primitive_symmetry(
        rotations, translations, kpoint, real=True
    )
    # Bilbao reports SM1, SM2 -- 2 PIRs of dim 2 each.
    assert len(physically_irreps) == 2

    table = get_cayley_table(pm_k_rotations)
    for pir in physically_irreps:
        assert pir.shape == (4, 2, 2)
        assert pir.dtype.kind == "f"
        assert is_representation(pir, table)
