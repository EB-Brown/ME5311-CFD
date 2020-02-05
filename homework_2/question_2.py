from typing import Callable, Dict, Iterable
from pathlib import Path

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

SIGMA = 11
RHO = 26
BETA = 8 / 3


def lorenz_system(point: np.array) -> np.array:
    """
    Calculate the lorenz-system vector.

    :param point: X, Y, Z coordinates
    :return: dx/dt, dy/dt, dz/dt
    """
    x, y, z = point
    return np.array([
        SIGMA * (y - x),
        x * (RHO - z) - y,
        x * y - BETA * z
    ])


def runge_kutta(func: Callable, initial: np.array, dt: float) -> np.array:
    """
    Vectorized Runge-Kutta numerical solver.

    :param func: First order differential equation
    :param initial: Initial condition f(0,0,0)
    :param dt: Time step
    :return: Numerical approximation of the next step
    """
    # Equation 2
    u_na = initial + dt * 8 / 15 * func(initial)

    # Equation 3
    u_na1 = u_na + dt * (-17 / 60 * func(initial) + 5 / 12 * func(u_na))

    # Equation 4
    return u_na1 + dt * (-5 / 12 * func(u_na) + 3 / 4 * func(u_na1))


def calc_array(num_points: int,
               initial: np.array,
               start_time: float,
               end_time: float) -> pd.DataFrame:
    """
    Calculates the array of approximated points for the ODE.

    :param num_points: Number of points along the interval
    :param initial: The initial condition in a `Point` object
    :param end_time: The end of the time interval
    :return: A list of points to be used as the numerical array
    """
    delta_time = (end_time - start_time) / num_points

    # initialize output with initial condition as first value
    point_array = {comp: [val] for comp, val in zip("xyz", initial)}
    point_array['t'] = [start_time]

    for i in range(1, 1 + num_points):
        # Store time stamp
        point_array['t'].append(start_time + i * delta_time)

        # Extract previous point from array
        previous_point = np.array([point_array[comp][-1] for comp in "xyz"])

        # Approximate next step
        xyz = runge_kutta(lorenz_system, previous_point, delta_time)

        # Store results
        for comp, val in zip('xyz', xyz):
            point_array[comp].append(val)

    return pd.DataFrame(point_array)


def vector_distances(a: pd.DataFrame,
                     b: pd.DataFrame) -> Dict[str, Iterable]:
    """
    Calculate the cartesian distances

    :param a: Control results
    :param b: Modified results
    :return: An array of the distance between the both set of results
    """
    assert ((a.t - b.t) == 0).all()
    return {'t': a.t, 'dist': np.sqrt(((a - b) ** 2).sum(axis=1))}


def plot_weather_prediction(results: pd.DataFrame, initial_condition: float):
    """
    Generates a plot of the weather prediction.

    :param results: A DataFrame of the weather prediction
    :param initial_condition:
        Initial condition used for analysis is injected into the plot title
    :return: Figure and axes subplot objects
    """
    # Initialize convergence plot
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.set_position([.15, .14, .575, .78])

    # Plot Components
    x_line, y_line, z_line = [
        ax.plot(experiment_control.t, vals)[0]
        for comp, vals in results.items()
        if comp != 't'
    ]

    ax.set_title(f"Weather Prediction: x0 = y0 = z0 = {initial_condition}")
    ax.set_xlabel("Time (Days)")
    ax.set_ylabel("Position (Miles)")

    plt.legend(
        [x_line, y_line, z_line],
        ['X Component', 'Y Component', 'Z Component'],
        loc='upper left',
        bbox_to_anchor=(1, 0.65),
    )
    return fig, ax


def lyapunov_func(t: float, lyapunov_exponent: float) -> float:
    """
    Standard function for the Lyapunov fit.

    :param t: time stamp
    :param t: Lyapunov exponent
    :return: Lyapunov function value
    """
    return np.exp(lyapunov_exponent * t)


if __name__ == "__main__":
    """Question 1"""
    initial_condition = np.array([1, 1, 1])  # x, y, z
    experiment_control = calc_array(
        150000, initial_condition, start_time=0, end_time=40
    )

    fig, ax = plot_weather_prediction(
        experiment_control, initial_condition[0]
    )

    output_dir = Path(__file__).parent / 'plots'
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / "Q2_1_weather.png")

    """Question 2, initial condition set to 1.001"""
    modified_initial_condition = initial_condition * 1.001
    experiment_modified = calc_array(
        150000, modified_initial_condition, start_time=0, end_time=40
    )

    fig, ax = plot_weather_prediction(
        experiment_modified, modified_initial_condition[0]
    )

    output_dir = Path(__file__).parent / 'plots'
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / "Q2_2_weather.png")

    """Plot Prediction Divergence"""
    distance = vector_distances(experiment_control, experiment_modified)
    (lyapunov_exp, *_), *_ = curve_fit(
        lyapunov_func, distance['t'], distance['dist']
    )
    fit = [lyapunov_func(t, lyapunov_exp) for t in distance['t']]

    # Initialize convergence plot
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.set_position([.15, .14, .575, .75])

    # Plot Components
    dist_line, *_ = ax.plot(distance['t'], distance['dist'])
    fit_line, *_ = ax.plot(distance['t'], fit)

    ax.set_title(
        f"Variation of Weather Prediction"
        f"\nLyapunov Exponent = {round(lyapunov_exp, 4)}"
    )
    ax.set_xlabel("Time (Days)")
    ax.set_ylabel("Variation Distance\n(sum([ai-bi]^2)^0.5 in Miles)")

    plt.legend(
        [dist_line, fit_line],
        ["Cartesian Distance", "Lyapunov Fitting"],
        loc='upper left',
        bbox_to_anchor=(1, 0.65),
    )

    output_dir = Path(__file__).parent / 'plots'
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / "Q2_2_variation.png")
