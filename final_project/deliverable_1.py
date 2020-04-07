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
    output_plots = Path(__file__).parent / 'plots'
    output_plots.mkdir(parents=True, exist_ok=True)
    for png in output_plots.rglob('*.png'):
        os.remove(png)

    # Analysis settings
    l_x = 2
    base_step = 2
    x_start, x_end = 0, 10
    y_start, y_end = 0, 10

    convergence_array = defaultdict(list)
    exponents = np.arange(2, 10.2, step=0.2)
    loops = len(exponents)
    step_size = exact_solution = mesh = None
    for progress, exponent in enumerate(exponents):
        step_size = 1 / base_step ** exponent

        # Initialize Mesh
        mesh = LinearMesh(
            x_start=x_start, x_end=x_end, x_step=step_size,
            y_start=y_start, y_end=y_end, y_step=step_size,
        )

        # Calculate numerical and exact solution
        second_derivative, exact_solution = d_x_y(mesh.x, mesh.y, l_x)
        numerical_solution = get_dct2_solution(
            mesh.x, mesh.y, second_derivative
        )

        # Store convergence_array information
        convergence_array['step_size'].append(step_size)
        convergence_array['inf_norm'].append(
            infinity_norm(numerical_solution, exact_solution)
        )

        # Make a contour plot of numerical solution
        if not (progress + 1) % 5:
            fig, ax = plot_contour(
                x_array=mesh.x,
                y_array=mesh.y,
                z_array=numerical_solution,
                title=f"Numerical Solution - dx=dy=1/{step_size}"
            )

            fig.savefig(
                output_plots / f"press_cont_{int(1 / step_size)}_points.png"
            )
        print(f"Completed progress: {round(100 * (1 + progress) / loops, 3)}%")

    fig, ax = plot_contour(
        x_array=mesh.x,
        y_array=mesh.y,
        z_array=exact_solution,
        title=f"Exact Solution - dx=dy=1/{step_size}"
    )

    fig.savefig(output_plots / f"exact_{round(1 / step_size)}_points.png")

    _, convergence_rate, convergence_fit = get_convergence_fit(
        convergence_array['step_size'], convergence_array['inf_norm'],
    )

    fig, ax = plot_convergence(
        x_array=convergence_array['step_size'],
        residuals=convergence_array['inf_norm'],
        convergence_fit=convergence_fit,
        title=f"Residual Convergence: Convergence Rate = {convergence_rate}",
        x_label='Step Size (dx=dy)',
        y_label='Residual Error',
    )

    fig.savefig(output_plots / "convergence_plot.png")
