"""ML_CVD_Model — default solver (FDM-based, replaces old erf solver).

This module provides backward-compatible API while using the improved
Crank-Nicolson finite difference solver under the hood.
"""

import warnings
from .fdm import ML_CVD_FDM

__all__ = ["ML_CVD_Model"]


class ML_CVD_Model:
    """Multilayer CVD diffusion model (FDM solver)."""

    def __init__(self):
        self.c1_res = None
        self.c2_res = None
        self._x_position = None
        self.D = []
        self.dt_used = []

    def __call__(self, filepath, dcal_type=0, dx=1e-3):
        solver = ML_CVD_FDM()
        solver(filepath, dcal_type, dx)
        self.c1_res = solver.c1_res
        self.c2_res = solver.c2_res
        self._x_position = solver.x_position
        self.D = solver.D
        self.dt_used = solver.dt_used

    @property
    def x_position(self):
        """Spatial coordinate grid (um)."""
        return self._x_position

    @x_position.setter
    def x_position(self, value):
        self._x_position = value

    @property
    def x_positon(self):
        warnings.warn(
            "x_positon is deprecated, use x_position instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self._x_position
