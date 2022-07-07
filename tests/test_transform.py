import numpy as np
from spglib import get_symmetry_dataset

from spgrep.transform import (
    get_primitive_transformation_matrix,
    transform_symmetry_and_kpoint,
    unique_primitive_symmetry,
)


def test_transformation(corundum_cell):
    # T in conventional basis
    kpoint = np.array([0, 1, 1 / 2])
    dataset = get_symmetry_dataset(corundum_cell)
    to_primitive = get_primitive_transformation_matrix(dataset["hall_number"])
    rotations, translations, primitive_kpoint = transform_symmetry_and_kpoint(
        to_primitive, dataset["rotations"], dataset["translations"], kpoint
    )
    primitive_rotations, primitive_translations, mapping = unique_primitive_symmetry(
        rotations, translations
    )

    # CDML coefficient
    assert np.allclose(primitive_kpoint, np.array([1 / 2, 1 / 2, -1 / 2]))

    # 12 unique symmetry operations with primitive cell in R-3c.
    assert primitive_rotations.shape == (12, 3, 3)
    assert primitive_translations.shape == (12, 3)
    assert set(mapping) == set(range(12))

    assert np.allclose([np.abs(np.linalg.det(r)) for r in primitive_rotations], 1)
