"""Shared utilities for diffusion_sicm solvers.

Common functions used by both the analytical erf solver and
the finite difference solver.
"""

import numpy as np

__all__ = ["interp_temperature", "compute_diffusion", "E_OVER_K"]

E_OVER_K = 11604.5  # e/k in K/eV


def interp_temperature(t: float, t_arr: np.ndarray, T_arr: np.ndarray) -> float:
    """Linear interpolation of temperature at time t."""
    if t <= t_arr[0]:
        return T_arr[0]
    if t >= t_arr[-1]:
        return T_arr[-1]
    idx = 0
    for i in range(len(t_arr) - 1):
        if t_arr[i] <= t < t_arr[i + 1]:
            idx = i
            break
    f = (t - t_arr[idx]) / (t_arr[idx + 1] - t_arr[idx])
    return T_arr[idx] + f * (T_arr[idx + 1] - T_arr[idx])


def compute_diffusion(
    temperature: float,
    d_coff: float,
    d_exp_c: float,
    d_temp_ref: float,
    dcal_type: int,
) -> float:
    """Compute diffusion coefficient D at a given temperature."""
    if dcal_type == 0:
        return d_coff * 1e8 * np.exp(-d_exp_c * E_OVER_K / temperature)
    else:
        return d_coff * np.exp(d_exp_c * (1.0 / temperature - 1.0 / d_temp_ref))
