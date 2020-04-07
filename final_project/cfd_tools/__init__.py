__all__ = [
    "get_dct2_solution",
    "p_norm",
    "infinity_norm",
    "get_convergence_fit",
    "plot_contour",
    "LinearMesh",
    "plot_convergence",
]

from .poisson_pressure import get_dct2_solution
from .convergence import p_norm, infinity_norm, get_convergence_fit
from .plots import plot_contour, plot_convergence
from .mesh import LinearMesh

__version__ = 1.0
