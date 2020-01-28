from typing import Union, List, Tuple, Any
from collections import namedtuple
from pathlib import Path
import math as m

from multiprocessing import Pool
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

REAL = Union[float, int]
INTERVAL_ARRAY = List[Tuple[float, float]]
Result = namedtuple('Result', ['numerical', 'samples'])


def func(x: REAL) -> float:
    """The original function f(x)"""
    return 0.5 * (abs(m.cos(4 * m.pi * x)) + x)


def analytic_solution() -> float:
    """Calculate the analytic solution for the integral of f(x) from 0 to 1"""
    vals = []
    for low, high in [(0, 0.125), (3 / 8, 5 / 8), (7 / 8, 1)]:
        # Integrate over the positive intervals of cos(4 * pi * x)
        vals.append(sin_func(high) - sin_func(low))

    for low, high in [(1 / 8, 3 / 8), (5 / 8, 7 / 8)]:
        # Integrate over the negative intervals of cos(4 * pi * x)
        vals.append(sin_func(low) - sin_func(high))

    return (sum(vals) + 0.5) / 2


def sin_func(x: REAL) -> float:
    """Helper function for calculating the sin term of the analytic solution"""
    return m.sin(4 * m.pi * x) / (4 * m.pi)


def monte_carlo_method(pairs: np.ndarray) -> float:
    """
    Calculate the integral using the Monte Carlo method.

    :param pairs:
        A 2-D array of random data points on the interval [0, 1] X [0, 1]
    :return: Numerical approximation for the integral
    """
    return sum(get_y_lt_func(pairs)) / len(pairs)


def get_y_lt_func(pairs: np.ndarray) -> np.array:
    """
    Helper function for determining which data points have y values are smaller
    than f(x).

    :param pairs: A 2-D array of random data points on the interval [0, 1] X [0, 1] :return:
        A boolean array indicating which points have y values less than f(x)
    """
    return np.array([y < func(x) for x, y in pairs])


def get_pairs(size: REAL) -> np.ndarray:
    """
    Randomly generating data points.

    :param size: Number of data points to generate.
    :return:
        A 2-D array of random data points on the interval [0, 1] X [0, 1]
    """
    return np.random.random((size, 2))


# Generate points along f(x)
X_ARRAY = np.arange(start=0, stop=1, step=0.001)
FUNC_POINTS = [func(x) for x in X_ARRAY]


def create_plot(pairs: np.ndarray) -> Tuple[Figure, Any]:
    """
    Generate a plot of the randomly generated data set. Points where y is
    greater than f(x) are marked by red dots. Points where y is less than f(x)
    are marked by green dots.

    :param pairs:
        A 2-D array of random data points on the interval [0, 1] X [0, 1]
    :return: Figure and axis objects
    """
    # Get a boolean array of points where y is less than f(x)
    idx = get_y_lt_func(pairs)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))  # Initialize subplot
    ax.set_position([.1, .125, .59, .8])  # Position axis

    # original function
    f_of_x = ax.plot(X_ARRAY, FUNC_POINTS, color='black')

    # "less than" points
    scat_lt = ax.plot(
        pairs[idx, 0], pairs[idx, 1], color='green', marker='o', ls=''
    )

    # "greater than" points
    scat_gt = ax.plot(
        pairs[~idx, 0], pairs[~idx, 1], color='red', marker='o', ls=''
    )

    # Plot labels
    ax.set_title(f"{len(pairs)} Samples")
    ax.set_xlabel("X Value")
    ax.set_ylabel("Y Value")

    # Insert plot legend
    plt.legend(
        [f_of_x[0], scat_lt[0], scat_gt[0]],  # plotted data
        ['f(x)=(|cos(4*pi*x)|+x)/2', 'y < f(x)', 'y >= f(x)'],  # labels
        loc='upper left',
        bbox_to_anchor=(1, 0.65),
    )

    return fig, ax


def execute_monte_carlo(size: REAL) -> float:
    return monte_carlo_method(get_pairs(size))


def process_sample_size(size: int) -> Result:
    n_pairs = get_pairs(size)
    fig, _ = create_plot(n_pairs)
    fig.savefig(output_dir / f"monte_carlo_plot_{size}_samples.png")

    return Result(monte_carlo_method(n_pairs), size)


if __name__ == "__main__":
    # Establish output directory
    base_dir = Path(__file__).parent
    output_dir = base_dir / 'plots'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine analytic solution
    analytic_answer = analytic_solution()
    print(
        f"The analytic solution for the function is {round(analytic_answer, 3)}"
    )

    # Process multiple sample sizes
    pool = Pool(4)
    results = pd.DataFrame(pool.map(process_sample_size, [100, 1000]))

    # Calculate errors
    results.loc[:, 'errors'] = abs(results.numerical - analytic_answer)

    # Initialize convergence plot
    fig, ax = plt.subplots()
    ax.set_position([.2, .125, .7, .8])

    # Plot residuals
    ax.loglog(results.samples, results.errors, color='black', marker='o')

    # Insert plot labels
    ax.set_title("Residual Convergence")
    ax.set_xlabel('Samples (N)')
    ax.set_ylabel('Residual Error |numerical - analytical|')

    # Save figure
    fig.savefig(output_dir / 'convergence_plot.png')

    error_1, error_2 = results.errors.iloc[0], results.errors.iloc[-1]
    sample_1, sample_2 = results.samples.iloc[0], results.samples.iloc[-1]

    convergence_rate = (error_2 - error_1) / (sample_2 - sample_1)

    print(f"In this trial, the convergence rate was {convergence_rate}")

    print(results)