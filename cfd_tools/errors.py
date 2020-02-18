def convergence_fit(num_points: int, a: float, k: float) -> float:
    """
    Calculate convergence points

    :param num_points: Number of samples used
    :param a: Constant multiplier
    :param k: Convergence rate
    :return: Lorenzo Gradient
    """
    return a * num_points ** (-k)