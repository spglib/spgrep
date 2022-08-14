import numpy as np
import pytest

from spgrep.group import get_cayley_table, get_little_group
from spgrep.irreps import enumerate_small_representations, enumerate_unitary_irreps
from spgrep.representation import is_representation
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
                rep_additional = np.zeros((2 * half_dim, 2 * half_dim), dtype=np.complex_)
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
