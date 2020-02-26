from typing import Callable, Tuple, Union
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DOMAIN = Tuple[float, float]
ARRAY = Union[np.array, pd.Series]

OUTPUT_DIR = Path(__file__).parent / 'plots'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def u_exact(x_array: ARRAY, t_step: float) -> np.array:
    """
    Exact solution to the convection equation given the initial condition.

    u_t + u_x = 0

    u(x, t) = f(x - t)
    u_t = f' * -1
    u_x = f' * 1
    u_t + u_x = f' * -1 + f' * 1 = 0

    Boundary condition:
    let, u(x, t) = cos (4*pi*[x - t])
    u(x, t=0) = cos (4*pi*[x - 0]) = cos(4*pi*x)

    :param x_array: Array of x values
    :param t_step: current timestamp
    :return: Flow velocity
    """
    return np.cos(4 * np.pi * (x_array - t_step))


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


def stencils(u_array: np.array,
             dx: float,
             one_step_coef: float,
             two_step_coef: float) -> np.array:
    """
    Approximates the first order derivative using spatial derivatives of two and
    one step away.

    This function applies periodic boundary conditions along the edge of the
    spatial domain.

    :param u_array: Array of u(x) values
    :param dx: time step
    :param one_step_coef: Coefficient for the one-step component
    :param two_step_coef: Coefficient for the one-step component
    :return: Approximated second derivative array
    """
    u_minus_2 = np.array([*u_array[-2:], *u_array[:-2]])
    u_plus_2 = np.array([*u_array[2:], *u_array[:2]])
    two_steps_out = two_step_coef * (u_plus_2 - u_minus_2)

    u_minus_1 = np.array([u_array[-1], *u_array[:-1]])
    u_plus_1 = np.array([*u_array[1:], u_array[0]])
    one_step_out = one_step_coef * (u_plus_1 - u_minus_1)

    return (two_steps_out + one_step_out) / dx


def run_simple_convection(x_domain: DOMAIN,
                          t_domain: DOMAIN,
                          dx: float,
                          dt: float,
                          one_step_coef: float,
                          two_step_coef: float) -> pd.DataFrame:
    """
    Run a simulation of the convection problem presented in question 2 of
    homework 4.

    :param x_domain: Lower and upper limits of the space domain
    :param t_domain: Lower and upper limits of the time domain
    :param dx: Step along the x component
    :param dt: Time step
    :param one_step_coef:
        Coefficient for the one-step component of the spatial derivative
    :param two_step_coef:
        Coefficient for the two-step component of the spatial derivative
    :return: A DataFrame of results from the simulation
    """
    x_start, x_end = x_domain
    x_array = np.arange(x_start, x_end + dx, dx)
    u_t0 = np.cos(4 * np.pi * x_array)  # u(t=0, x)

    t_start, t_end = t_domain
    t_array = np.arange(t_start, t_end + dt, dt)[1:]  # Skip initial condition

    # New function with fixed dx and coefficients
    def du_dx(previous_step: ARRAY) -> ARRAY:
        """A helper function for flipping the sign on the first derivative"""
        return -1 * stencils(
            previous_step,
            dx=dx,
            one_step_coef=one_step_coef,
            two_step_coef=two_step_coef,
        )

    output = {0: u_t0}  # Initial condition
    old_t = 0
    for t in t_array:
        # Last calculated step
        last_time_step = output[old_t]

        # Calculate the next step
        # Note: pu/pt = - pu/px = spatial derivative
        current_time_step = runge_kutta(du_dx, last_time_step, dt)

        # Drop NaNs and save to output
        output[t] = current_time_step
        old_t = t

    # Format data types
    output = pd.DataFrame(output, index=x_array)

    # Return only values within the
    return output


def triple_comparison_plot(x_array: ARRAY,
                           u_exact: ARRAY,
                           u_q1a: ARRAY,
                           u_q1b: ARRAY,
                           title: str):
    """
    Generate the plot for question 2a.

    This will compare the exact solution to the two finite difference methods
    presented in question 1.

    :param x_array: Spatial Array
    :param u_exact:
        Exact solution for the convection equation at a given time step
    :param u_q1a: Numerical solution from equation 1a
    :param u_q1b: Numerical solution from equation 1b
    :param title: Plot title
    :return: figure and axis objects from matplotlib.pyplot.subplots
    """
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.set_position([.15, .14, .575, .78])
    exact, *_ = ax.plot(x_array, u_exact, color="black")
    q1a_numeric, *_ = ax.plot(x_array, u_q1a, color='red', ls='--')
    q1b_numeric, *_ = ax.plot(x_array, u_q1b, color='green', ls='dotted')
    ax.set_title(title)
    ax.set_xlabel("X Position")
    ax.set_ylabel(" Flow Velocity")
    ax.set_ylim(-1.2, 1.2)
    plt.legend(
        [exact, q1a_numeric, q1b_numeric],
        ['Exact Solution', "Equation 1", "Equation 2"],
        loc='upper left',
        bbox_to_anchor=(1, 0.65),
    )
    return fig, ax


def comparison_plot(x_array: ARRAY,
                    u_exact: ARRAY,
                    u_numeric: ARRAY,
                    title: str):
    """
    Helper function for generating plots comparing the exact vs numerical
    solution.

    :param x_array: X array
    :param u_exact: Array of the exact solution
    :param u_numeric: Array of the numeric solution
    :param title: Plot title
    :return: Figure and axis object
    """
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.set_position([.15, .14, .575, .78])
    ex, *_ = ax.plot(x_array, u_exact, color='black')
    num, *_ = ax.plot(x_array, u_numeric, color='red', ls='--')
    ax.set_title(title)
    ax.set_xlabel("X Position")
    ax.set_ylabel("Velocity")
    plt.legend(
        [ex, num],
        ['Exact Solution', 'Equation 2'],
        loc='upper left',
        bbox_to_anchor=(1, 0.65),
    )
    return fig, ax


def contour_plot(df: pd.DataFrame, title: str):
    """
    Generate a contour plot of the solution for the convection equation

    :param df: Resulting values for the convection equation
    :param title: Title of the plot
    :return: Figure and axis object
    """
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
    """Question 2 part 1"""
    dx = 0.05
    dt = dx / 10
    x_domain = 0, 1
    t_domain = 0, 1
    x_array = np.arange(x_domain[0], x_domain[1] + dx, dx)

    exact_solution = pd.DataFrame({
        round(t, 5): u_exact(x_array, t)
        for t in np.arange(t_domain[0], t_domain[1] + dt, dt)
    })
    exact_solution.index = x_array

    solution_1a = run_simple_convection(
        x_domain=x_domain,
        t_domain=t_domain,
        dx=dx,
        dt=dt,
        one_step_coef=2 / 3,
        two_step_coef=- 1 / 12,
    )

    solution_1b = run_simple_convection(
        x_domain=x_domain,
        t_domain=t_domain,
        dx=dx,
        dt=dt,
        one_step_coef=1,
        two_step_coef=- 1 / 4,
    )

    question_2a = {
        "Exact Solution": exact_solution,
        "Equation 1": solution_1a,
        "Equation 2": solution_1b,
    }

    question_2_a_output = OUTPUT_DIR / 'question_2_1'
    question_2_a_output.mkdir(parents=True, exist_ok=True)

    for t in np.arange(0, 1.01, 0.1):
        t_1a = np.where(solution_1a.columns >= t)[0][0]
        t_1b = np.where(solution_1b.columns >= t)[0][0]
        t_u_exact = np.where(exact_solution.columns >= t)[0][0]

        timestamp_1a = solution_1a.columns[t_1a]
        timestamp_1b = solution_1b.columns[t_1b]
        timestamp_exact = exact_solution.columns[t_u_exact]

        fig, ax = triple_comparison_plot(
            x_array=solution_1a.index,
            u_exact=exact_solution[timestamp_exact],
            u_q1a=solution_1a[timestamp_1a],
            u_q1b=solution_1b[timestamp_1b],
            title=f"t = {timestamp_1b}"
        )
        ts = str(round(timestamp_1b, 1)).replace(".", "_")
        fig.savefig(question_2_a_output / f"Solution Plot - t_{ts}.png")
        plt.close(fig)

    for method, results in question_2a.items():
        fig, ax = contour_plot(results, title=f"Velocity Propagation: {method}")
        fig.savefig(question_2_a_output / f"{method} Contour Plot.png")
        plt.close(fig)

    """Question 2 part 2"""
    question_2_b_output = OUTPUT_DIR / 'question_2_2'
    question_2_b_output.mkdir(parents=True, exist_ok=True)
    for dx_val in [0.05, 0.025, 0.0123, 0.00625]:
        dt_val = dx_val / 10

        x2_array = np.arange(x_domain[0], x_domain[1] + dx_val, dx_val)
        exact = pd.DataFrame({
            round(t, 5): u_exact(x2_array, t)
            for t in np.arange(t_domain[0], t_domain[1] + dt_val, dt_val)
        })
        exact.index = x2_array

        numeric = run_simple_convection(
            x_domain=x_domain,
            t_domain=t_domain,
            dx=dx_val,
            dt=dt_val,
            one_step_coef=1,
            two_step_coef=-1 / 4,
        )

        t_1b = np.where(numeric.columns >= t_domain[1])[0][0]
        timestamp_1b = numeric.columns[t_1b]

        t_u_exact = np.where(exact.columns >= t_domain[1])[0][0]
        timestamp_exact = exact.columns[t_u_exact]

        fig, ax = comparison_plot(
            x2_array,
            exact[timestamp_exact],
            numeric[timestamp_1b],
            f"Numerical vs Exact Solution: dx = {dx_val}",
        )
        fig.savefig(question_2_b_output / f"Comparison Plot {dx_val}.png")
        plt.close(fig)
