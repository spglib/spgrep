import numpy as np

from spgrep.utils import (
    contain_space,
    is_integer_array,
    ndarray2d_to_integer_tuple,
    nroot,
)


def test_is_integer_array():
    assert is_integer_array(np.array([1.0, 2.0, 4.0]))
    assert not is_integer_array(np.array([[1.0, 2.0, 4.0], [1.5, 2.5, 4.5]]))


def test_ndarray2d_to_integer_tuple():
    actual = ndarray2d_to_integer_tuple(np.array([[-1.0, 0.0], [0.0, 1.0]]))
    expect = ((-1, 0), (0, 1))
    assert actual == expect


def test_nroot():
    actual = nroot(-1j, 2)
    expect = np.exp(-1j * np.pi * 1 / 4)
    assert np.allclose(actual, expect)


def test_contain_space():
    # Strictly contained
    assert contain_space(
        np.array([[1, 0, 0], [0, 1, 0]]),
        np.array([[1, 1, 0]]),
    )

    # Equal case
    assert contain_space(
        np.array([[1, 0, 0], [0, 1, 0]]),
        np.array([[1, 1, 0], [1, -1, 0]]),
    )

    # False cases
    assert not contain_space(
        np.array([[1, 0, 0], [0, 1, 0]]),
        np.array([[1, 0, 1]]),
    )
    assert not contain_space(
        np.array([[1, 0, 0], [0, 1, 0]]),
        np.array([[1, 0, 0], [0, 0, 1]]),
    )
    assert not contain_space(
        np.array([[1, 0, 0], [0, 1, 0]]),
        np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    )
