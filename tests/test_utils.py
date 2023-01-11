import numpy as np
import pytest

from spgrep.utils import (
    grassmann_distance,
    is_integer_array,
    is_prime,
    mode_dot,
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


@pytest.mark.parametrize(
    "n,expect",
    [
        (42, False),
        (43, True),
    ],
)
def test_is_prime(n, expect):
    assert is_prime(n) == expect


def test_nroot():
    actual = nroot(-1j, 2)
    expect = np.exp(-1j * np.pi * 1 / 4)
    assert np.allclose(actual, expect)


def test_grassmann_distance():
    assert (
        grassmann_distance(
            np.array(
                [
                    [1, 0, 0],
                    [0, 1j, 0],
                ]
            ),
            np.array(
                [
                    [0, 1, 0],
                    [0, 0, 1],
                ]
            ),
        )
        > 0
    )

    assert np.isclose(
        grassmann_distance(
            np.array(
                [
                    [2, 0, 0],
                    [0, 4j, 0],
                ]
            ),
            np.array(
                [
                    [1, 1j, 0],
                    [8, 0, 0],
                ]
            ),
        ),
        0,
        atol=1e-5,
    )

    assert np.isclose(
        grassmann_distance(
            np.array(
                [
                    [1, 1j, 0],
                ]
            ),
            np.array(
                [
                    [0, 1, 0],
                    [1, 0, 0],
                ]
            ),
        ),
        0,
        atol=1e-5,
    )


def test_mode_dot(rng):
    # Test 3-mode product
    coeffs = rng.random((2, 2, 2))
    m1 = rng.random((2, 3))
    m2 = rng.random((2, 5))
    m3 = rng.random((2, 7))

    expect = np.einsum("ijk,ia,jb,kc", coeffs, m1, m2, m3)
    actual = mode_dot(coeffs, [m1, m2, m3])
    assert actual.shape == (3, 5, 7)
    assert np.allclose(actual, expect)
