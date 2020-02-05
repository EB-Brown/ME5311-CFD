from typing import Callable, List
from collections import namedtuple, defaultdict
from pathlib import Path

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd

Point = namedtuple("Point", ['x', 'y'])


def runge_kutta(func: Callable, initial: float, dt: float) -> float:
    """
    Runge-Kutta numerical solver.

    :param func: First order differential equation
    :param initial: Initial condition
    :param dt: Time step
    :return: Numerical approximation of the next step
    """
    # Equation 2
    u_na = initial + dt * 8 / 15 * func(initial)

    # Equation 3
    u_na1 = u_na + dt * (-17 / 60 * func(initial) + 5 / 12 * func(u_na))

    # Equation 4
    return u_na1 + dt * (-5 / 12 * func(u_na) + 3 / 4 * func(u_na1))


def ode(u: float) -> float:
    """
    The given ODE (Equation 5).

    :param u: Input value for u(t)
    :return: Value of u'(t)
    """
    return u ** 2


def exact_func(t: float) -> float:
    """
    The exact function u(t) to be compared against the numerical calculations.

    :param t: Time value
    :return: u(t)
    """
    return 1 / (1 - t)


def calc_array(num_points: int, initial: Point, end_time: float) -> List[Point]:
    """
    Calculates the array of approximated points for the ODE.

    :param num_points: Number of points along the interval
    :param initial: The initial condition in a `Point` object
    :param end_time: The end of the time interval
    :return: A list of points to be used as the numerical array
    """
    delta_time = (end_time - initial.x) / num_points
    point_array = [initial]
    for i in range(1, 1 + num_points):
        x_val = initial.x + i * delta_time
        point_array.append(
            Point(
                x=x_val,
                y=runge_kutta(  # Approximate next value
                    func=ode,  # Equation 5
                    initial=point_array[-1][1],  # Previously calculated value
                    dt=delta_time
                ),
            )
        )

    return point_array


def exact_array(num_points: int,
                start_time: float,
                end_time: float) -> List[Point]:
    """
    Calculate an array of exact values.

    :param num_points: Number of points along the interval
    :param start_time: The start of the time interval
    :param end_time: The end of the time interval
    :return: A list of points to be used as the numerical array
    """
    delta_time = (end_time - start_time) / num_points
    point_array = []
    for i in range(1 + num_points):
        x_val = start_time + i * delta_time
        point_array.append(
            Point(x=x_val, y=exact_func(x_val))
        )
    return point_array


def p_norm(numerical: List[Point], exact: List[Point], order: int = 2) -> float:
    """
     Calculate a p_norm of the error.

    :param numerical: An array of numerically approximated values
    :param exact: An array of the exact values
    :param order: Order of the p-norm summation
    :return: Norm Error
    """
    dt = (exact[-1][0] - exact[0][0]) / len(exact)
    errors = []
    for n, ((_, aprx), (_, exct)) in enumerate(zip(numerical, exact)):
        errors.append(abs(aprx ** n - exct ** n) ** order)

    return (dt * sum(errors)) ** (1 / order)


def convergence_fit(num_points: int, a: float, k: float) -> float:
    """
    Calculate convergence points

    :param num_points: Number of samples used
    :param a: Constant multiplier
    :param k: Convergence rate
    :return: Lorenzo Gradient
    """
    return a * num_points ** (-k)


if __name__ == '__main__':
    initial_condition = Point(2, -1)
    p_norms = defaultdict(list)
    for x in range(1, 101):
        number_of_points = 150 * x

        numerical_solver = calc_array(
            number_of_points, initial_condition, end_time=3
        )
        exact = exact_array(number_of_points, start_time=2, end_time=3)

        p_norms['number_of_points'].append(number_of_points)
        p_norms['p-norms'].append(
            p_norm(numerical=numerical_solver, exact=exact)
        )

    (a, convergence_rate), *_ = curve_fit(
        convergence_fit, p_norms['number_of_points'], p_norms['p-norms']
    )

    fit = [
        convergence_fit(n, a, convergence_rate)
        for n in p_norms['number_of_points']
    ]

    # Initialize convergence plot
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.set_position([.15, .14, .575, .78])

    # Plot residuals
    experiment, *_ = ax.loglog(
        p_norms['number_of_points'],
        p_norms['p-norms'],
        color='red',
        marker='x',
        ls='',
    )

    # Plot residual fitting
    residual_fit, *_ = ax.loglog(p_norms['number_of_points'], fit, color='k')

    # Insert plot labels
    convergence_rate = round(convergence_rate, 3)
    ax.set_title(f"Residual Convergence: Convergence Rate = {convergence_rate}")
    ax.set_xlabel('Samples (log scale)\n(N)')
    ax.set_ylabel('2nd Order Error Norm (log scale)')

    # insert plot legend
    plt.legend(
        [residual_fit, experiment],
        ["Convergence Fit", "Experiment Residuals"],
        loc='upper left',
        bbox_to_anchor=(1, 0.65),
    )

    # Save figure
    output_dir = Path(__file__).parent / 'plots'
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / 'Q1_convergence_plot.png')

    """Solver Plot"""
    # Initialize solver plot
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.set_position([.15, .14, .575, .78])

    # Reformat arrays
    numerical_solver = pd.DataFrame(numerical_solver)
    exact = pd.DataFrame(exact)

    # Plot arrays
    exact_line, *_ = ax.plot(exact.x, exact.y, color='black')
    numerical, *_ = ax.plot(
        numerical_solver.x, numerical_solver.y, color='red', ls='--'
    )

    # Insert plot labels
    ax.set_title(f"Numerical vs Exact Solution")
    ax.set_xlabel('X Value')
    ax.set_ylabel('Y Value')

    # insert plot legend
    plt.legend(
        [numerical, exact_line],
        ["Numerical Solution", "Exact Solution"],
        loc='upper left',
        bbox_to_anchor=(1, 0.65),
    )

    # Save figure
    output_dir = Path(__file__).parent / 'plots'
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / 'Q1_numerical_vs_exact.png')

    fig, ax = plt.subplots()
