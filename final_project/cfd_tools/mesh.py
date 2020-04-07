from collections import defaultdict

import numpy as np


class LinearMesh:

    def __init__(self, **kwargs):
        parameters = defaultdict(dict)
        for key, val in kwargs.items():
            *parm, arg = key.split("_")
            parameters["_".join(parm)][arg] = val

        for parm, domain in parameters.items():
            start = domain.get('start', 0)
            end = domain.get('end', start + 1)
            step = domain.get('step', (end - start) / 10)
            extended_end = end + step

            array = np.arange(start, extended_end, step)
            array = array[array <= end]
            setattr(self, parm, array)
