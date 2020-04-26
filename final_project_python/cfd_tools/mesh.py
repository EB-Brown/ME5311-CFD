__all__ = ['LinearMesh']

from collections import defaultdict

import numpy as np


class LinearMesh:

    def __init__(self, **kwargs):
        """
        This class is used to generate a N-dimension hyper mesh.

        Each axis can be dynamically named. The axis must include a `_start`,
        `_end`, and `_step` arguement.

        Example usage:

            mesh = LinearMesh(
                x_start=1, x_end=2, x_step=0.1,
                y_start=3, y_end=4, y_step=0.2,
            )

        """
        # Parse inputs
        parameters = defaultdict(dict)
        for key, val in kwargs.items():
            *parm, arg = key.split("_")
            parameters["_".join(parm)][arg] = val

        # Generate mesh arrays
        for parm, domain in parameters.items():
            start = domain.get('start', 0)
            end = domain.get('end', start + 1)
            step = domain.get('step', 0)
            extended_end = end + step

            if step != 0:  # Generate an array
                array = np.arange(start, extended_end, step)
                array = array[array <= end]
            else:  # Generate an array with a single value
                array = np.array([start])

            setattr(self, parm, array)
            setattr(self, f"{parm}_step", step)
