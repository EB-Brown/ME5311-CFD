from typing import Callable, Tuple
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DOMAIN = Tuple[float, float]

OUTPUT_DIR = Path(__file__).parent / 'plot'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def runge_kutta(func: Callable, initial: np.array, dt: float) -> np.array:
    """
    Runge-Kutta numerical solver.

    :param func: First order differential equation
    :param initial: Initial condition
    :param dt: Time step
    :return: Numerical approximation of the next step
    """
    u_na = initial + dt * 8 / 15 * func(initial)
    u_na1 = u_na + dt * (-17 / 60 * func(initial) + 5 / 12 * func(u_na))
    return u_na1 + dt * (-5 / 12 * func(u_na) + 3 / 4 * func(u_na1))


def three_point_stencil(u_array: np.array, dx: float) -> np.array:
    """
    Approximates the second order derivative using the three point stencil
    centered finite difference

    :param u_array: Array of u(x) values
    :param dx: time step
    :return: Approximated second derivative array
    """
    return np.array([0, *((u_array[:-2] - u_array[2:]) / (2 * dx)), 0])


def run_simple_convection(x_domain: DOMAIN,
                          t_domain: DOMAIN,
                          dx: float,
                          dt: float) -> pd.DataFrame:
    """
    Run a simulation of the convection problem presented in question 2 of
    homework 4.

    :param x_domain: Lower and upper limits of the space domain
    :param t_domain: Lower and upper limits of the time domain
    :param dx: Step along the x component
    :param dt: Time step
    :return: A DataFrame of results from the simulation
    """
    t_start, t_end = t_domain
    t_array = np.arange(t_start, t_end + dt, dt)[1:]  # Skip initial condition
    u_x0 = np.sin(t_array * 4)  # u(t, x=0)

    # Each iteration drops the last data point in the x domain
    # To compensate, x_array extends to x_end + dx * len(t_array)
    x_start, x_end = x_domain
    x_array = np.arange(x_start, x_end + dx * len(t_array), dx)
    u_t0 = np.zeros(len(x_array))  # u(t=0, x)

    du_dx = partial(three_point_stencil, dx=dx)  # New function with fixed dx

    output = {0: u_t0}  # Initial condition
    for t, u_at_x_equals_0 in zip(t_array, u_x0):
        # Last calculated step
        last_time_step = output[round(t - dt, 5)]

        # Calculate the next step
        # Note: pu/pt = - pu/px = three_point_stencil
        current_time_step = runge_kutta(du_dx, last_time_step, dt)
        # Insert boundary at x = 0 to first value of array
        current_time_step.itemset(0, u_at_x_equals_0)

        # Drop NaNs and save to output
        output[round(t, 5)] = current_time_step[~np.isnan(current_time_step)]

        if t >= t_end:
            break

    # Format data types
    output = pd.DataFrame(output)

    # Set index to x values
    output.index = x_array[:len(output)]

    # Return only values within the
    return output.loc[output.index <= t_end]


def line_plot(x, y, title):
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.set_position([.15, .14, .575, .78])
    ax.plot(x, y)
    ax.set_title(title)
    ax.set_xlabel("X Position")
    ax.set_ylabel("Velocity")
    return fig, ax


def contour_plot(df: pd.DataFrame, title: str):
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.set_position([.15, .14, .575, .78])
    for plot_grain in reversed(range(1000)):
        try:  # Plot with max possible contour lines
            contour = ax.contour(
                list(df.columns),
                df.index,
                df.values,
                plot_grain,
                cmap='RdGy'
            )
            break  # It worked!
        except RuntimeError:  # Did not work
            pass

    ax.set_title(title)
    ax.set_ylabel("X Position")
    ax.set_xlabel("Time")
    plt.colorbar(contour)

    return fig, ax


if __name__ == '__main__':
    dx = 0.01
    dt = dx * 3 ** 0.5
    solution = run_simple_convection(
        x_domain=(0, 1),
        t_domain=(0, 2),
        dx=dx,
        dt=dt,
    )

    stable_output = OUTPUT_DIR / 'stable_solution'
    stable_output.mkdir(parents=True, exist_ok=True)

    for t in np.arange(0, 2.01, 0.25):
        t = np.where(solution.columns >= t)[0][0]
        timestamp = solution.columns[t]

        fig, ax = line_plot(
            x=solution.index,
            y=solution.iloc[:, t],
            title=f"t = {timestamp}"
        )
        ts = str(timestamp).replace(".", "_")
        fig.savefig(stable_output / f"Solution Plot - t_{ts}.png")

    fig, ax = line_plot(
        x=solution.columns,
        y=solution.iloc[0],
        title="Initial Condition: x = 0"
    )
    fig.savefig(stable_output / "Initial Condition - x_0.png")

    fig, ax = contour_plot(
        solution,
        title=f"Velocity Propagation: dx = {dx}, dt = {round(dt, 5)}"
    )
    fig.savefig(stable_output / "Contour Plot.png")

    """Unstable Analysis"""

    unstable = run_simple_convection(
        x_domain=(0, 1),
        t_domain=(0, 2),
        dx=dx,
        dt=dt * 1.1,
    )

    unstable_output = OUTPUT_DIR / 'unstable_solution'
    unstable_output.mkdir(parents=True, exist_ok=True)

    for t in np.arange(0, 2, 0.25):
        t = np.where(unstable.columns >= t)[0][0]
        timestamp = solution.columns[t]

        fig, ax = line_plot(
            x=unstable.index,
            y=unstable.iloc[:, t],
            title=f"t = {unstable.columns[t]}"
        )

        ts = str(timestamp).replace(".", "_")
        fig.savefig(unstable_output / f"Solution Plot - t_{ts}.png")

    fig, ax = contour_plot(
        unstable,
        title=f"Velocity Propagation: dx = {dx}, dt = {round(dt * 1.1, 5)}"
    )
    fig.savefig(unstable_output / "Contour Plot.png")

    print()
