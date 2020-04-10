__all__ = ['p_norm', 'infinity_norm', 'get_convergence_fit']

from typing import Union
from collections import namedtuple

import numpy as np
from scipy.optimize import curve_fit

converg_fit = namedtuple("converg_fit", 'a rate fit')


def p_norm(numerical: np.ndarray,
           exact: np.ndarray,
           step_size: float,
           order: int = 2) -> float:
    """
    Calculate a p_norm of the error for a 1-D problem.

    :param numerical: An array of numerically approximated values
    :param exact: An array of the exact values
    :param step_size: Step size for each solution
    :param order: Order of the p-norm summation
    :return: Norm Error
    """
    return (step_size * (abs(numerical - exact) ** order).sum()) ** (1 / order)


def infinity_norm(numerical: np.ndarray,
                  exact: np.ndarray) -> float:
    """
    Calculate the maximum error between a numerical solution and the exact
    solution.

    :param numerical: An array of numerically approximated values
    :param exact: An array of the exact values
    :return: Infinity norm of the numerical solution
    """
    return abs(numerical - exact).max()


def _step_convergence(x_array: Union[int, np.ndarray],
                      a: float,
                      k: float) -> float:
    """
    A helper function for calculating the convergence points given the step size
    used for each analytic.

    :param x_array: Array of step sizes
    :param a: Constant multiplier
    :param k: Convergence rate
    :return: Lorenzo Gradient
    """
    return a * x_array ** k


def _point_convergence(x_array: Union[int, np.ndarray],
                       a: float,
                       k: float) -> float:
    """
    A helper function for calculating the convergence points given the number of
    points used for each analytic.

    :param x_array: Number of samples used
    :param a: Constant multiplier
    :param k: Convergence rate
    :return: Lorenzo Gradient
    """
    return _step_convergence(x_array, a, -1 * k)


def get_convergence_fit(x_array: np.ndarray,
                        y_array: np.ndarray,
                        step_fit: bool = True) -> converg_fit:
    """
    Calculates the convergence fit properties and the fit array.

    :param x_array:
    :param y_array:
    :param step_fit:
    :return:
        `namedtuple` containing the convergence-fit multiplier,
        convergence_rate, and the fitted values for the given `x_array`
    """
    fitter = _step_convergence if step_fit else _point_convergence
    (a, convergence_rate), *_ = curve_fit(fitter, x_array, y_array)

    return converg_fit(
        a, convergence_rate, fitter(x_array, a=a, k=convergence_rate)
    )
