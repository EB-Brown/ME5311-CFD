__all__ = ['FluidProfile', 'Field']

from abc import ABC
from typing import Tuple

import numpy as np

from cfd_tools import LinearMesh, get_dct2_solution

VELOCITY_TUPLES = Tuple[np.ndarray, np.ndarray]


class Field:
    values: np.ndarray
    mesh: LinearMesh
    _partial_mesh: LinearMesh = None
    _dy: np.ndarray = None
    _dx: np.ndarray = None

    def __init__(self,
                 mesh: np.ndarray,
                 values: LinearMesh):
        """"""
        self.values = values
        self.mesh = mesh

    @property
    def dx(self) -> np.ndarray:
        if self._dx is None:
            self._dx = \
                (self.values[:, 1:] - self.values[:, :-1]) / self.mesh.x_step
        return self._dx

    @property
    def dy(self) -> np.ndarray:
        if self._dy is None:
            self._dy = \
                (self.values[1:] - self.values[:-1]) / self.mesh.y_step
        return self._dy

    @property
    def partial_mesh(self) -> LinearMesh:
        if self._partial_mesh is None:
            self._partial_mesh = LinearMesh(
                x_start=self.mesh.x[0] + self.mesh.x_step / 2,
                x_end=self.mesh.x[-1] - self.mesh.x_step / 2,
                x_step=self.mesh.x_step,
                y_start=self.mesh.y[0] + self.mesh.y_step / 2,
                y_end=self.mesh.y[-1] - self.mesh.y_step / 2,
                y_step=self.mesh.y_step,
            )
        return self._partial_mesh


class FluidProfile(ABC):
    _x_velocity: Field
    _y_velocity: Field
    __pressure: Field = None

    def __init__(self, x_velocity: Field, y_velocity: Field):
        """"""
        self._x_velocity: Field = x_velocity
        self._y_velocity: Field = y_velocity

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
    def _pressure(self) -> Field:
        if self.__pressure is None:
            laplacian_pressure = self.velocity_divergence
            self.__pressure = self._poisson_solver(laplacian_pressure)
        return self.__pressure

    @property
    def pressure(self) -> np.ndarray:
        return self._pressure.values

    @property
    def x_velocity(self) -> np.ndarray:
        """check divergence, correct, then send"""
        return self._x_velocity.values

    @property
    def y_velocity(self) -> np.ndarray:
        """check divergence, correct, then send"""
        return self._y_velocity.values

    def remove_divergence(self):
        # Get velocity profiles
        x_velocity = self._x_velocity.values
        y_velocity = self._y_velocity.values

        # Remove divergence
        x_velocity[:, 1:-1] -= self._pressure.dx
        y_velocity[1:-1] -= self._pressure.dy

        # Assert boundaries
        x_velocity, y_velocity = \
            self._assert_boundary_conditions(x_velocity, y_velocity)

        # Update velocity profiles
        self._x_velocity = Field(mesh=self._x_velocity.mesh, values=x_velocity)
        self._y_velocity = Field(mesh=self._y_velocity.mesh, values=y_velocity)

    def _assert_boundary_conditions(self,
                                    x_velocity: np.ndarray,
                                    y_velocity: np.ndarray) -> VELOCITY_TUPLES:
        """"""
        y_velocity[:2] = 0
        y_velocity[-2:] = 0

        x_velocity[:, :2] = 0
        x_velocity[:, -2:] = 0

        return x_velocity, y_velocity

    def _poisson_solver(self, laplacian_pressure: np.ndarray) -> Field:
        """
        Integrate the poisson problem for the pressure gradient.

        The returned results will be shifted a half step in the X and Y axis.

        :param laplacian_pressure:
        :return:
        """
        pressure = get_dct2_solution(
            x_array=self.pressure_x_domain,
            y_array=self.pressure_y_domain,
            z_array=laplacian_pressure,
        )

        return Field(mesh=self.pressure_mesh, values=pressure)
