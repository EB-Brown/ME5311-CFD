from typing import List
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
from scipy.optimize import curve_fit

POISSON_AMPLITUDE = -4 * np.pi ** 2
COS_FREQUENCY = 2 * np.pi


def exact_solution(x_array: np.ndarray) -> np.ndarray:
    return np.cos(COS_FREQUENCY * x_array)


def d_x(x_array: np.ndarray) -> np.ndarray:
    return POISSON_AMPLITUDE * exact_solution(x_array)


def get_dct_solution(x_array: np.ndarray, dx: float) -> np.ndarray:
    x_len = len(x_array)
    d = d_x(x_array)

    cos_trans = fftpack.dct(d, norm='ortho')  # ortho = MATLAB func
    denominator = (2 * np.cos(np.pi * (np.arange(x_len)) / x_len) - 2)
    p_k = (dx ** 2) * cos_trans / denominator

    # Handle first wave number where zero is in the denominator
    idx = np.isinf(p_k)
    p_k[idx] = 0

    return fftpack.idct(p_k, norm='ortho')


def p_norm(numerical: np.ndarray,
           exact: np.ndarray,
           dx: float,
           order: int = 2) -> float:
    """
     Calculate a p_norm of the error.

    :param numerical: An array of numerically approximated values
    :param exact: An array of the exact values
    :param dx: Step size
    :param order: Order of the p-norm summation
    :return: Norm Error
    """
    return (dx * sum(abs(numerical - exact) ** order)) ** (1 / order)


def convergence_fit(num_points: int, a: float, k: float) -> float:
    """
    Calculate convergence points

    :param num_points: Number of samples used
    :param a: Constant multiplier
    :param k: Convergence rate
    :return: Lorenzo Gradient
    """
    return a * num_points ** k


if __name__ == "__main__":
    output_dir = Path(__file__).parent / 'plots/question_2'
    output_dir.mkdir(parents=True, exist_ok=True)

    p_norms = {}
    for magnitude in range(1, 5):
        for i in range(1, 11):
            step = 1 / (i * 10 ** magnitude)
            x_vals = np.arange(0, 2, step)
            num_pressure = get_dct_solution(x_vals, step)
            exact_pressure = exact_solution(x_vals)
            error_norm = p_norm(num_pressure, exact_pressure, step)
            p_norms[step] = error_norm

        fig, ax = plt.subplots(figsize=(8, 4.8))
        ax.set_position([.15, .14, .575, .78])

        exact_line, *_ = ax.plot(x_vals, exact_pressure, color='black')
        num_line, *_ = ax.plot(x_vals, num_pressure, color='red', ls='--')

        ax.set_title(f"Numerical vs Exact Solution")
        ax.set_xlabel('X Position')
        ax.set_ylabel('Pressure')

        plt.legend(
            [num_line, exact_line],
            ['DCT Approximation', 'Analytical Solution'],
            loc='upper left',
            bbox_to_anchor=(1, 0.65),
        )
        fig.savefig(output_dir / f'dct_{round(1 / step)}_points.png')

    p_norms = pd.Series(p_norms).dropna().sort_index()
    (a, convergence_rate), *_ = curve_fit(
        convergence_fit, p_norms.index, p_norms
    )
    fit = [convergence_fit(n, a, convergence_rate) for n in p_norms.index]

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.set_position([.15, .14, .575, .78])

    experiment, *_ = ax.loglog(
        p_norms.index, p_norms, color='red', marker='x', ls='',
    )
    residual_fit, *_ = ax.loglog(p_norms.index, fit, color='k')

    convergence_rate = round(convergence_rate, 3)
    ax.set_title(f"Residual Convergence: Convergence Rate = {convergence_rate}")
    ax.set_xlabel('Step (Delta X)\nlog scale')
    ax.set_ylabel('2nd Order Error Norm (log scale)')

    plt.legend(
        [residual_fit, experiment],
        ["Convergence Fit", "Experiment Residuals"],
        loc='upper left',
        bbox_to_anchor=(1, 0.65),
    )

    fig.savefig(output_dir / 'convergence_plot.png')

    print()
