from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


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


def u_func(x: np.array, t: np.array):
    theta = x - 2 * t
    theta = np.where((0 < theta) & (theta <= 1), theta, 2 - theta)
    return np.where((0 < theta) & (theta <= 2), theta, 0)


if __name__ == "__main__":
    x_array = np.arange(-0.5, 4.5, 0.01)
    t_array = np.arange(0, 2, 0.01)
    solution = pd.DataFrame({t: u_func(x_array, t) for t in t_array})
    solution.index = x_array

    OUTPUT_DIR = Path(__file__).parent / 'plot'
    question_1 = OUTPUT_DIR / 'question_1'
    question_1.mkdir(parents=True, exist_ok=True)

    for t in np.arange(0, 2, 0.25):
        t = np.where(solution.columns >= t)[0][0]
        fig, ax = line_plot(
            x=solution.index,
            y=solution.iloc[:, t],
            title=f"t = {solution.columns[t]}"
        )
        fig.savefig(question_1 / f"Line Plot - t_{t}.png")

    fig, ax = contour_plot(
        solution,
        title=f"Velocity Propagation: dx = {0.01}, dt = {0.01}"
    )
    fig.savefig(question_1 / "Contour Plot.png")