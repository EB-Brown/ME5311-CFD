from typing import Dict, Callable, Tuple, List
from collections import defaultdict
from pathlib import Path
from functools import partial

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy

OUTPUT_DIR = Path(__file__).parent / 'plot'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def u_x(x: np.array) -> np.array:
    """
    The original function

    :param x: An array of X values
    :return: An array of U values
    """
    return x ** 2 - 7 / 8 * x ** 4 + 7 / 16 * x ** 6 - 3 / 16 * x ** 7 \
           + 3 / 128 * x ** 8


def u_x_ddot_exact(x: np.array) -> np.array:
    """
    Calculate the exact solution to the second derivative

    :param x: An array of x values
    :return: An array of the second derivatives
    """
    return 2 - 21 / 2 * x ** 2 + 105 / 8 * x ** 4 - 126 / 16 * x ** 5 \
           + 21 / 16 * x ** 6


def three_point_stencil(u_array: np.array, dx: float) -> np.array:
    """
    Approximates the second order derivative using the three point stencil
    centered finite difference

    :param u_array: Array of u(x) values
    :param dx: time step
    :return: Approximated second derivative array
    """
    u_i_minus_1 = u_array[:-2]
    u_i = u_array[1:-1]
    u_i_plus_1 = u_array[2:]
    return (u_i_minus_1 - 2 * u_i + u_i_plus_1) / (dx ** 2)


def convergence_fit(num_points: int, a: float, k: float) -> float:
    """
    Calculate convergence points

    :param num_points: Number of samples used
    :param a: Constant multiplier
    :param k: Convergence rate
    """
    return a * num_points ** k


def get_convergence(dt_array: np.array,
                    residuals: np.array) -> Tuple[float, Callable]:
    """
    Calculate the convergence rate and return a curve fit function

    :param dt_array: An array indicating the number of points in each test
    :param residuals: P Norm error for each test
    """
    (a, convergence_rate), *_ = curve_fit(
        convergence_fit, dt_array, residuals,
    )
    return convergence_rate, partial(convergence_fit, a=a, k=convergence_rate)


def execute_study(dx: float, interval_end: float) -> Dict[str, np.array]:
    """
    Generates the exact and approximations for the second derivative.

    :param dx: time step
    :param interval_end:
        End value of the interval. (2 for the homework assignment)
    """
    # Set up
    number_of_points = round(interval_end / dx)
    x = np.arange(0, interval_end, dx)

    # Exact solutions
    u_exact = u_x(x)
    u_ddot = u_x_ddot_exact(x)

    # Three point stencil
    three_point = three_point_stencil(u_exact, dx)

    # Spectral approximation
    fft = scipy.fft.fft(u_exact)
    kx = np.zeros(number_of_points)
    wave_length = 2 * np.pi / interval_end
    kx[:number_of_points // 2 + 1] = \
        np.arange(number_of_points // 2 + 1) * wave_length
    kx[number_of_points // 2 + 1:] = \
        np.arange(-1 * number_of_points // 2 + 1, 0) * wave_length
    u_ddot_fft = scipy.fft.ifft(-(kx ** 2) * fft).real

    return {'x': x, 'exact': u_ddot, 'stencil': three_point, 'fft': u_ddot_fft}


def p_norm(x_array: np.array,
           numerical: np.array,
           exact: np.array,
           order: int = 1) -> float:
    """
     Calculate a p_norm of the error.

    :param x_array: An array of x values
    :param numerical: An array of numerically approximated values
    :param exact: An array of the exact values
    :param order: Order of the p-norm summation
    :return: Norm Error
    """
    dt = (x_array[-1] - x_array[0]) / len(x_array)
    errors = []
    for n, (aprx, exct) in enumerate(zip(numerical, exact)):
        errors.append(abs(aprx ** n - exct ** n) ** order)

    return (dt * sum(errors)) ** (1 / order)


def generate_test_plots(tests: List[Dict[str, np.array]]):
    """
    Generate plots of each approximation fitting to the exact second derivative
    for a specific method.

    :param tests: Trial results for a specifc method
    """
    for test in tests:
        generate_test_plot(test['x'], test['exact'], test['fft'], 'fft')
        generate_test_plot(test['x'], test['exact'], test['stencil'], 'stencil')


TITLES = {
    'fft': 'Spectral Aprx.',
    'stencil': '3 Point Stencil'
}


def generate_test_plot(x: np.array, u: np.array, fit: np.array, method: str):
    """
    Generates a plot of the fittings.

    :param x: X array
    :param u: Exact solution array
    :param fit: Approximated array
    :param method: A string indicating the method type
    """
    # Initialize solver plot
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.set_position([.15, .14, .575, .78])

    exact, *_ = ax.plot(x, u, color='black')
    fit, *_ = ax.plot(
        x[1:-1] if len(x) != len(fit) else x,
        fit,
        color='red',
        ls='--'
    )

    dt = round((x[-1] - x[0]) / len(x), 4)
    title = TITLES[method]
    ax.set_title(title + f": dx={dt}")
    ax.set_ylabel('d2u/dx2')
    ax.set_xlabel('X position')
    ax.set_ylim(-2.5, 2.5)

    plt.legend(
        [exact, fit],
        ['Exact Solution', title],
        loc='upper left',
        bbox_to_anchor=(1, 0.65),
    )

    fig.savefig(OUTPUT_DIR / (title + f"dx_{dt}".replace(".", "_") + ".png"))


def generate_convergence_plots(norm: pd.DataFrame, method: str):
    """
    Generate a plot showing the convergence of the errors.

    :param norm: Residual norms
    :param method: Method used to approximate the second derivative
    """
    p_norm_reversed = list(reversed(norm['p_norm']))
    con_rate, fit = get_convergence(norm['dt'], p_norm_reversed)
    residual_fit = fit(norm['dt'])  # norm['num_points'])

    # Initialize solver plot
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.set_position([.15, .14, .575, .78])

    residuals, *_ = ax.loglog(
        norm['dt'], p_norm_reversed, color='black', marker='o', ls=''
    )
    fit, *_ = ax.loglog(norm['dt'], residual_fit, color='red')

    title = TITLES[method]
    ax.set_title(title + f": Convergence Rate = {con_rate}")
    ax.set_ylabel('First Order P Norm')
    ax.set_xlabel('dt')

    plt.legend(
        [residuals, fit],
        ['Residual Error', 'Convergence Fit'],
        loc='upper left',
        bbox_to_anchor=(1, 0.65),
    )

    fig.savefig(OUTPUT_DIR / (title + "_convergence.png"))


if __name__ == '__main__':
    # Run studies
    experiments = [
        execute_study(dx=1 / x, interval_end=2)
        for x in [8, 16, 32, 64, 128, 256]
    ]

    # Plot fittings
    generate_test_plots(experiments)

    # Calculate errors
    norms = defaultdict(list)
    for t in experiments:
        dt = (t['x'][-1] - t['x'][0]) / len(t['x'])
        num_points = len(t['x'])
        norms['fft'].append({
            'num_points': num_points,
            'dt': dt,
            'p_norm': p_norm(t['x'], t['fft'], t['exact'])
        })
        norms['stencil'].append({
            'num_points': num_points,
            'dt': dt,
            'p_norm': p_norm(t['x'][1:-1], t['stencil'], t['exact'][1:-1])
        })

    # Plot convergence
    # for method, norms in norms.items():
    #     generate_convergence_plots(pd.DataFrame(norms), method)
