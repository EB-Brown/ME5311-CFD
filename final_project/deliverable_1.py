from typing import Tuple

import numpy as np
from scipy import fftpack

SIZE_CACHE = {}


def _get_size(array: np.ndarray, var: str) -> Tuple[float, float]:
    """A helper function for determining if the domain is a constant distance"""
    val_id = id(array)
    try:  # Pull from the cache to avoid calculating a double diff
        return SIZE_CACHE[val_id]
    except KeyError:
        pass

    # Check is steps are consistent
    diff = np.diff(array)
    if not np.isclose(np.diff(diff), 0).all():
        raise ValueError(f"'{var}' does not have consistent time steps")

    output = SIZE_CACHE[val_id] = len(array), diff.mean()
    return output


KMOD_CACHE = {}


def _get_k_mod(x_len: float, y_len: float, dx: float, dy: float) -> np.ndarray:
    val_id = (id(arg) for arg in [x_len, y_len, dx, dy])
    try:  # Pull from the cache to avoid calculating the k mod array again
        return KMOD_CACHE[val_id]
    except KeyError:
        pass

    return (
        (2 * np.cos(np.pi * np.arange(x_len) / x_len) - 2) / dx ** 2 +
        (2 * np.cos(np.pi * np.arange(y_len) / y_len) - 2) / dy ** 2
    )


def get_dct2_solution(x_array: np.ndarray,
                      y_array: np.ndarray,
                      z_array: np.ndarray) -> np.ndarray:
    """


    :param x_array:
        Array of X-values. The array must have constant step distances.
    :param y_array:
        Array of Y-values. The array must have constant step distances.
    :param z_array:
        A 2-D array with X values as the columns and the Y values as the rows.
    :return: A double integral of the inputs array.
    """
    x_len, dx = _get_size(x_array, "x_array")
    y_len, dy = _get_size(y_array, "y_array")
    cos_trans = fftpack.fft2(z_array).real
    denominator = _get_k_mod(x_len, y_len, dx, dy)
    p_k = cos_trans / denominator

    # Handle first wave number where zero is in the denominator
    idx = np.isinf(p_k)
    p_k[idx] = 0

    return fftpack.ifft2(p_k).real / 2


def d_x_y(x_array: np.ndarray,
          y_array: np.ndarray,
          l_x: float) -> Tuple[np.ndarray, np.ndarray]:
    """

    :param x_array:
    :param y_array:
    :param l_x:
    :return: Second derivative and the exact solution
    """
    dd_array = []
    z_array = []
    for y in y_array:
        dd_y_calcs = []
        y_calcs = []
        for x in x_array:
            dd_y_calcs.append(
                -16 * np.pi ** 2 * (1 + 4 / l_x ** 2)
                * np.cos(8 * np.pi * x / l_x)
                * np.cos(4 * np.pi * y)
            )
            y_calcs.append(
                np.cos(8 * np.pi * x / l_x) * np.cos(4 * np.pi * y)
            )
        dd_array.append(dd_y_calcs)
        z_array.append(y_calcs)

    return np.array(dd_array), np.array(z_array)
