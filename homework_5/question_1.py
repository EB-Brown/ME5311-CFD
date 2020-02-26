from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

OUTPUT_DIR = Path("plots")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def q1a_k_mod(k_delta_x: np.array) -> np.array:
    """


    :param k_delta_x:
    :return:
    """
    return (8 * np.sin(k_delta_x) - np.sin(2 * k_delta_x)) / 6


def q1b_k_mod(k_delta_x: np.array) -> np.array:
    """


    :param k_delta_x:
    :return:
    """
    return (2 * np.sin(k_delta_x) - np.sin(2 * k_delta_x)) / 2


def calc_phase_speed(k_delta_x: np.array, k_mod: np.array) -> np.array:
    """


    :param k_delta_x:
    :param k_mod:
    :return:
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


def generate_plots(k_delta_x, k_mod, phase_speed, question_number: int):
    """


    :param k_delta_x:
    :param k_mod:
    :param phase_speed:
    :param question_number:
    question_number
    """
    question_directory = OUTPUT_DIR / f"question_1_{question_number}"
    question_directory.mkdir(parents=True, exist_ok=True)

    title = f"Q1.{question_number} - "
    fig, _ = plot_mod_wavenum(
        k_delta_x,
        k_mod,
        title=title + "Modified Wavenumber"
    )
    fig.savefig(question_directory / "modified_wavenumber.png")

    fig, _ = plot_phase_speed(
        k_delta_x,
        phase_speed,
        title=title + "Numerical Phase Speed"
    )
    fig.savefig(question_directory / "numerical_phase_speed.png")


if __name__ == "__main__":

    k_delta_x = np.arange(0, np.pi, 0.01)
    for q_num, func in enumerate((q1a_k_mod, q1b_k_mod)):
        values = {"k_delta_x": k_delta_x}
        values['k_mod'] = k_mod = func(k_delta_x)
        values['phase_speed'] = phase_speed = calc_phase_speed(
            k_delta_x, k_mod
        )
        generate_plots(**values, question_number=q_num + 1)
