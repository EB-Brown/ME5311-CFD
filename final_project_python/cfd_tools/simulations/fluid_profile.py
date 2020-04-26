__all__ = ['FluidProfile']

from abc import ABC
from typing import Tuple

import numpy as np

from .fields import Velocity as Vel, Pressure
from cfd_tools import LinearMesh, get_dct2_solution

VELOCITY_TUPLES = Tuple[np.ndarray, np.ndarray]


class FluidProfile(ABC):
    _x_velocity: Vel
    _y_velocity: Vel
    __pressure: Pressure = None

    def __init__(self, x_velocity: Vel, y_velocity: Vel):
        """"""
        self._x_velocity = x_velocity
        self._y_velocity = y_velocity

    @property
    def velocity_divergence(self) -> np.ndarray:
        return self._x_velocity.dx + self._y_velocity.dy

    @property
    def pressure_x_domain(self) -> np.ndarray:
        return self._x_velocity.partial_mesh.x

    @property
    def pressure_y_domain(self) -> np.ndarray:
        return self._y_velocity.partial_mesh.y

    @property
    def pressure_mesh(self) -> LinearMesh:
        return self._y_velocity.partial_mesh

    @property
    def _pressure(self) -> Pressure:
        if self.__pressure is None:
            laplacian_pressure = self.velocity_divergence
            self.__pressure = self._poisson_solver(laplacian_pressure)
        return self.__pressure

    @property
    def pressure(self) -> np.ndarray:
        return self._pressure.values

    @property
    def x_velocity(self) -> np.ndarray:
        """check velocity_divergence.m, correct, then send"""
        return self._x_velocity.values

    @property
    def y_velocity(self) -> np.ndarray:
        """check velocity_divergence.m, correct, then send"""
        return self._y_velocity.values

    def remove_divergence(self):
        # Get velocity profiles
        x_velocity = self._x_velocity.values
        y_velocity = self._y_velocity.values

        # Remove velocity_divergence.m
        x_velocity[1:-2, 1:-1] -= self._pressure.dx
        y_velocity[1:-1, 1:-2] -= self._pressure.dy

        # Assert boundaries.m
        x_velocity, y_velocity = \
            self._assert_boundary_conditions(x_velocity, y_velocity)

        # Update velocity profiles
        self._x_velocity = Vel(mesh=self._x_velocity.mesh, values=x_velocity)
        self._y_velocity = Vel(mesh=self._y_velocity.mesh, values=y_velocity)

    @staticmethod
    def _assert_boundary_conditions(x_velocity: np.ndarray,
                                    y_velocity: np.ndarray) -> VELOCITY_TUPLES:
        """"""
        # Bottom wall
        y_velocity[0] = 0
        # x_velocity[0] = 0

        # Top wall
        y_velocity[-2] = 0
        # x_velocity[-2] = 0

        # Left wall
        # y_velocity[:, 0] = 0
        x_velocity[:, 0] = 0

        # Right wall
        # y_velocity[:, -1] = 0
        x_velocity[:, -1] = 0

        return x_velocity, y_velocity

    def _poisson_solver(self, laplacian_pressure: np.ndarray) -> Pressure:
        """
        Integrate the poisson problem for the pressure gradient.

        The returned results will be shifted a half step in the X and Y axis.

        :param laplacian_pressure:
        :return:
        """
        pressure = get_dct2_solution(
            x_array=self.pressure_x_domain[:-1],
            y_array=self.pressure_y_domain[:-1],
            z_array=laplacian_pressure,
        )

        return Pressure(mesh=self.pressure_mesh, values=pressure)
