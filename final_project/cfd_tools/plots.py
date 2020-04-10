__all__ = ['plot_contour', 'plot_convergence']

import numpy as np
import matplotlib.pyplot as plt


def plot_contour(*,
                 x_array: np.ndarray,
                 y_array: np.ndarray,
                 z_array: np.ndarray,
                 title: str = None,
                 x_label: str = None,
                 y_label: str = None) -> plt.subplot:
    """
    Generate a contour plot of the solution for the convection equation

    :param x_array: First axis of the solution
    :param y_array: Second axis of the solution
    :param z_array: Third axis of the solution
    :param title: Title of the plot
    :param x_label: X-axis label
    :param y_label: Y-axis label
    :return: Figure and axis object
    """
    fig, ax = _get_plot()

    x_len = len(x_array)
    y_len = len(y_array)
    if x_len == 1:
        x_array = np.array([*x_array, *(x_array + 1)])
    if y_len == 1:
        y_array = np.array([*y_array, *(y_array + 1)])
    if 1 in (y_len, x_len):
        z_array = np.array([*z_array, *z_array])

    for plot_grain in reversed(range(1000)):
        try:  # Plot with max possible contour lines
            contour = ax.contour(
                x_array, y_array, z_array, plot_grain, cmap='RdGy'
            )
            break  # It worked!
        except RuntimeError:  # Did not work
            pass

    # Plot Legend
    plt.colorbar(contour)

    return fig, _write_labels(ax, title, x_label, y_label)


def plot_convergence(*,
                     x_array: np.ndarray,
                     residuals: np.ndarray,
                     convergence_fit: np.ndarray,
                     title: str = None,
                     x_label: str = None,
                     y_label: str = None) -> plt.subplot:
    """
    Generate a log-log plot of the residual error and rate of convergence.


    :param x_array: An array of the number of points or step sizes
    :param residuals: The solution's raw error
    :param convergence_fit: A logarithmic fit of the solution's error
    :param title: Title of the plot
    :param x_label: X-axis label
    :param y_label: Y-axis label
    :return: Figure and axis object
    """
    fig, ax = _get_plot()

    resid_points, *_ = ax.loglog(
        x_array, residuals, color='red', marker='o', ls='',
    )
    fit, *_ = ax.loglog(x_array, convergence_fit, color='k')
    plt.legend(
        [resid_points, fit],
        ["Residuals", "Convergence Fit"],
        loc='upper left',
        bbox_to_anchor=(1, 0.65),
    )
    return fig, _write_labels(ax, title, x_label, y_label)


def _get_plot() -> plt.subplot:
    """Helper function for initializing plot objects"""
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.set_position([.15, .14, .575, .78])
    return fig, ax


def _write_labels(axis,
                  title: str = None,
                  x_label: str = None,
                  y_label: str = None) -> plt.subplot:
    """Helper function for writing labels on the plots"""
    axis.set_title(title)
    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    return axis
