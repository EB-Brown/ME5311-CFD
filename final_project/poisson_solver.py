import os
from typing import Tuple
from pathlib import Path
from collections import defaultdict

import numpy as np

from cfd_tools import (
    LinearMesh,
    get_dct2_solution,
    infinity_norm,
    get_convergence_fit,
    plot_contour,
    plot_convergence,
)


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

    try:
        y_array[0]
    except TypeError:
        y_array = [y_array]

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
    output_plots = Path(__file__).parent / 'plots/poisson_solver'
    output_plots.mkdir(parents=True, exist_ok=True)
    for png in output_plots.rglob('*.png'):
        os.remove(png)

    # Analysis settings
    base_step = 2
    x_start, x_end = 0, 1.5
    y_start, y_end = 0, 2
    x_len = x_end - x_start
    y_len = y_end - y_start

    convergence_array = defaultdict(list)
    exponents = np.arange(2, 10, step=0.2)
    loops = len(exponents)
    step_size = exact_solution = mesh = numerical_solution = num_points = None
    for progress, exponent in enumerate(exponents):
        num_points = base_step ** exponent
        x_step_size = x_len / num_points
        y_step_size = y_len / num_points
        step_size = 1 / num_points

        # Initialize Mesh
        mesh = LinearMesh(
            x_start=x_start, x_end=x_end, x_step=x_step_size,
            y_start=y_start, y_end=y_end, y_step=y_step_size,
        )

        # Calculate numerical and exact solution
        second_derivative, exact_solution = d_x_y(mesh.x, mesh.y, l_x=x_len)
        numerical_solution = get_dct2_solution(
            mesh.x, mesh.y, second_derivative
        )

        # Store convergence_array information
        convergence_array['step_size'].append(step_size)
        convergence_array['inf_norm'].append(
            infinity_norm(numerical_solution, exact_solution)
        )

        # Make a contour plot of numerical solution
        if not (progress + 1) % 3 or (progress + 1) == loops:
            fig, ax = plot_contour(
                x_array=mesh.x,
                y_array=mesh.y,
                z_array=numerical_solution,
                title=f"Numerical Solution - dx=dy=1/{num_points}"
            )

            fig.savefig(
                output_plots / f"press_cont_{int(num_points)}_points.png"
            )
        print(f"Completed progress: {round(100 * (1 + progress) / loops, 3)}%")

    num_points = int(num_points)
    fig, ax = plot_contour(
        x_array=mesh.x,
        y_array=mesh.y,
        z_array=exact_solution,
        title=f"Exact Solution - Number of Points = {num_points}"
    )

    fig.savefig(output_plots / f"exact_{num_points}_points.png")

    perc_err = \
        100 * abs(exact_solution - numerical_solution) / exact_solution.max()
    resid_title = "Percent Residual Error Solution - Number of Points = {}"
    fig, ax = plot_contour(
        x_array=mesh.x,
        y_array=mesh.y,
        z_array=perc_err,
        title=resid_title.format(num_points)
    )

    fig.savefig(output_plots / f"residuals_{num_points}_points.png")

    _, convergence_rate, convergence_fit = get_convergence_fit(
        convergence_array['step_size'], convergence_array['inf_norm'],
    )
    convergence_rate = round(convergence_rate, 3)

    fig, ax = plot_convergence(
        x_array=convergence_array['step_size'],
        residuals=convergence_array['inf_norm'],
        convergence_fit=convergence_fit,
        title=f"Residual Convergence: Convergence Rate = {convergence_rate}",
        x_label='Step Size (dx=dy)',
        y_label='Residual Error',
    )

    fig.savefig(output_plots / "convergence_plot.png")
