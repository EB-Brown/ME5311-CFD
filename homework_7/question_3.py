from typing import Mapping
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags

POISSON_AMPLITUDE = -4 * np.pi ** 2
COS_FREQUENCY = 2 * np.pi


def exact_solution(x_array: np.ndarray) -> np.ndarray:
    return np.cos(COS_FREQUENCY * x_array)


def d_x(x_array: np.ndarray) -> np.ndarray:
    return POISSON_AMPLITUDE * exact_solution(x_array)


def jacoby_iteration(b_i: pd.Series,
                     scheme: np.ndarray,
                     initial_guess: pd.Series,
                     relax: float = 1) -> np.array:
    l = len(scheme)
    assert l % 2  # Handles only odd length schemes
    diagonal = scheme[l // 2]
    scheme[l // 2] = 0  # Avoid summing diagonal term

    # Rolling calculation
    a_x = initial_guess.rolling(l).apply(lambda x: sum(scheme * x))

    # Rolling calculation assigns val to shifted position
    a_x = pd.concat([a_x.iloc[1:], a_x.iloc[0:1]])
    a_x.index = b_i.index

    return relax / diagonal * (b_i - a_x) + (1 - relax) * initial_guess


def get_jacoby_solution(b_i: pd.Series,
                        scheme: np.ndarray,
                        position_conditions: Mapping,
                        initial_guess: pd.Series,
                        relax: float = 1,
                        num_of_iterations: int = 100):
    y = initial_guess.copy()
    b_mag = sum(b_i.iloc[1:-1] ** 2) ** 0.5
    b_slice = b_i.iloc[1:-1]

    residuals = []
    for i in range(num_of_iterations):
        y = jacoby_iteration(b_i, scheme.copy(), y, relax)

        # Insert boundary conditions
        for x, val in position_conditions.items():
            y.loc[x] = val

        # Calculate residuals
        y_m1 = y.iloc[:-2].reset_index(drop=True)
        y_m1.index = b_slice.index
        y_0 = y.iloc[1:-1].reset_index(drop=True)
        y_0.index = b_slice.index
        y_p1 = y.iloc[2:].reset_index(drop=True)
        y_p1.index = b_slice.index
        ddy = scheme[0] * y_m1 + scheme[1] * y_0 + scheme[2] * y_p1

        r = ((sum((ddy - b_slice) ** 2) ** 0.5) / b_mag)

        residuals.append({'iteration': i + 1, 'residual': r})

        if r > 100:
            break

    return y, pd.DataFrame(residuals)


def question_3(start: float,
               end: float,
               dx: float,
               relax: float = 1,
               num_of_iterations: int = 100):
    x_array = np.arange(start, end + dx, dx)
    b_i = pd.Series(d_x(x_array))
    b_i.index = x_array

    scheme = np.array([1, -2, 1]) / dx ** 2
    boundary_conditions = {start: 1, end: 1}

    initial_guess = pd.Series(np.zeros(len(x_array)))
    initial_guess.index = b_i.index
    for x, val in boundary_conditions.items():
        initial_guess.loc[x] = val

    y, residuals = get_jacoby_solution(
        b_i,
        scheme,
        boundary_conditions,
        initial_guess,
        relax,
        num_of_iterations
    )
    return x_array, y, residuals


if __name__ == "__main__":
    output_dir = Path(__file__).parent / 'plots/question_3'
    output_dir.mkdir(parents=True, exist_ok=True)

    dx = 0.02
    size = 100
    scheme = np.array([1, -2, 1]) / dx
    A = diags(scheme, offsets=[-1, 0, 1], shape=(size, size)).toarray()

    _, eig_vec = np.linalg.eig(A)
    p = np.dot(np.dot(np.linalg.inv(eig_vec), A), eig_vec)
    p_inv = np.linalg.inv(p)

    jacoby_matrix = np.identity(size) - np.dot(p_inv, A)
    rho = np.linalg.eigvals(jacoby_matrix).max()
    optimal = 1 + (rho / (1 + (1 + rho * 2) ** 0.5) ** 2) ** 2

    exact = None
    residual_df = {}
    for r in [0.5, 0.75, 1, optimal]:
        x, y, resid = question_3(
            start=0, end=2, dx=dx, relax=r, num_of_iterations=2000
        )

        residual_df[r] = resid

        if exact is None:
            exact = exact_solution(x)

        fig, ax = plt.subplots(figsize=(8, 4.8))
        ax.set_position([.15, .14, .575, .78])
        ex, *_ = ax.plot(x, exact, color='black')
        num, *_ = ax.plot(x, y, color='red', ls='--')
        ax.set_title(f"Jacobi Fit: Omega = {round(r, 3)}")
        ax.set_xlabel('X Position')
        ax.set_ylabel('Pressure')
        plt.legend(
            [ex, num],
            ['Exact Solution', 'Jacoby Solution'],
            loc='upper left',
            bbox_to_anchor=(1, 0.65),
        )
        if r == optimal:
            r = 'optimal'
        else:
            r = f'{round(r, 3)}'.replace(".", "_")
        fig.savefig(output_dir / f'jacobi_fit_omega_{r}.png')

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.set_position([.15, .14, .575, .78])

    lines = []
    for r, df in residual_df.items():
        lin, *_ = ax.loglog(df.iteration, df.residual)
        if r == optimal:
            title = 'Optimal Relaxation'
        else:
            title = f"Relaxation = {round(r, 3)}"
        lines.append((lin, title))

    ax.set_title(f"Relaxation Convergence")
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Residuals')
    plt.legend(
        [l for l, _ in lines],
        [t for _, t in lines],
        loc='upper left',
        bbox_to_anchor=(1, 0.65),
    )
    fig.savefig(output_dir / "relaxation_residuals.png")
    plt.show()

    print()
