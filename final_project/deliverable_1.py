from typing import Tuple

import numpy as np

from cfd_tools import get_dct2_solution, infinity_norm


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


if __name__ == "__main__":
    # TODO
    #  define base step size
    #  For loop: iterable factor to decrease step size
    #    create x_array and y_array
    #    calculate second derivative and exact solution
    #    feed to dct2_solver
    #    Get infinity norm
    #    plot solutions
    #  get convergence fit
    #  generate convergence plot

    print()
