import numpy as np
from spglib import get_symmetry_dataset

from spgrep.transform import (
    get_primitive_transformation_matrix,
    transform_symmetry_and_kpoint,
)


def test_transformation(corundum_cell):
    # T in conventional basis
    kpoint = np.array([0, 1, 1 / 2])
    dataset = get_symmetry_dataset(corundum_cell)
    to_primitive = get_primitive_transformation_matrix(dataset["hall_number"])
    primitive_rotations, primitive_translations, primitive_kpoint = transform_symmetry_and_kpoint(
        to_primitive, dataset["rotations"], dataset["translations"], kpoint
    )

    # CDML coefficient
    assert np.allclose(primitive_kpoint, np.array([1 / 2, 1 / 2, -1 / 2]))
