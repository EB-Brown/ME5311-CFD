from typing import Dict, Callable
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUTPUT_DIR = Path("plots")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def hw_5_1a_k_mod(k_delta_x: np.array) -> np.array:
    """
    Modified wave number for the standard fourth-order explicit scheme of
    Homework # 5.

    :param k_delta_x: array representing k delta x
    :return: Modified wave number * delta x
    """
    return (8 * np.sin(k_delta_x) - np.sin(2 * k_delta_x)) / 6


A = 17 / 12
B = 101 / 150
C = 1 / 100
ALPHA = 1 / 2
BETA = 1 / 20


def hw6_k_mod(k_delta_x: np.array) -> np.array:
    """
    Modified wave number for the tenth-order implicit scheme in homework 6.

    :param k_delta_x: array representing k delta x
    :return: Modified wave number * delta x
    """
    return (
               2 * C * np.sin(3 * k_delta_x)
               + 3 * B * np.sin(2 * k_delta_x)
               + 6 * A * np.sin(k_delta_x)
           ) / (
               (
                   2 * BETA * np.cos(2 * k_delta_x)
                   + 2 * ALPHA * np.cos(k_delta_x)
                   + 1
               )
               * 6
           )


def calc_phase_speed(k_delta_x: np.array, k_mod: np.array) -> np.array:
    """
    Calculate the numerical-phase velocity (Implicit filtering property).


    :param k_delta_x: array representing k delta x
    :param k_mod: Modified wave number * delta x
    :return: The numerical-phase velocity (Implicit filtering property)
    """
    return k_mod / k_delta_x


def plot_mod_wavenum(k_delta_x: np.array, k_mod_delta_x: np.array, title: str):
    """


    :param k_delta_x:
    :param k_mod_delta_x:
    :param title:
    :return:
    """
    fig, ax = plt.subplots()
    ax.plot(k_delta_x, k_mod_delta_x)
    ax.set_title(title)
    ax.set_xlabel(r"$k \Delta x$")
    ax.set_ylabel(r"$k^{*} \Delta x$")
    return fig, ax


def plot_phase_speed(k_delta_x: np.array, phase_speed: np.array, title: str):
    """

    :param k_delta_x:
    :param phase_speed:
    :param title:
    :return:
    """
    fig, ax = plt.subplots()
    ax.plot(k_delta_x, phase_speed)
    ax.set_title(title)
    ax.set_xlabel(r"$k \Delta x$")
    ax.set_ylabel(r"$a^{*}/a$")
    return fig, ax


def generate_plots(k_delta_x,
                   k_mod,
                   phase_speed,
                   scheme_name: str,
                   ppw: int):
    """


    :param k_delta_x:
    :param k_mod:
    :param phase_speed:
    :param scheme_name:
    :param ppw: Points per Wavelength
    """
    question_directory = OUTPUT_DIR / f"{scheme_name.replace(' ', '_').lower()}"
    question_directory.mkdir(parents=True, exist_ok=True)

    title = f"{scheme_name} - "
    fig, _ = plot_mod_wavenum(
        k_delta_x,
        k_mod,
        title=title + "Modified Wave Number"
    )
    fig.savefig(question_directory / "modified_wavenumber.png")

    fig, _ = plot_phase_speed(
        k_delta_x,
        phase_speed,
        title=title + f"Phase Velocity, PPW = {int(ppw)}"
    )
    fig.savefig(question_directory / "phase_velocity.png")

    if "Ten" in scheme_name:
        fig, _ = plot_phase_speed(
            k_delta_x,
            phase_speed,
            title=title + f"Implicit Filter, PPW = {int(ppw)}"
        )
        fig.savefig(question_directory / "implicit_filter.png")


if __name__ == "__main__":

    k_delta_x = np.arange(0, np.pi, 0.001)
    f_order = "Fourth-Order Explicit Scheme"
    t_order = "Tenth-Order Implicit Scheme"
    max_error = 0.01  # 1% error
    scheme: Dict[str, Callable] = {f_order: hw_5_1a_k_mod, t_order: hw6_k_mod}

    values = {}
    for func_name, func in scheme.items():
        values[func_name] = {"k_delta_x": k_delta_x}
        values[func_name]['k_mod'] = k_mod = func(k_delta_x)
        values[func_name]['phase_speed'] = phase_speed = \
            calc_phase_speed(k_delta_x, k_mod)

        idx = phase_speed >= 1 - max_error
        ppw = 2 * np.pi // k_mod[idx].max()
        generate_plots(**values[func_name], scheme_name=func_name, ppw=ppw)
        values[func_name] = pd.DataFrame(values[func_name])
