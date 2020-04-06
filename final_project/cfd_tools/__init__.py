__all__ = [
    "get_dct2_solution",
    "p_norm",
    "infinity_norm",
    "get_convergence_fit",
]

from .poisson_pressure import get_dct2_solution
from .convergence import p_norm, infinity_norm, get_convergence_fit

__version__ = 1.0
