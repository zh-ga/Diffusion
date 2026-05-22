"""Finite Difference Method solver for multilayer CVD diffusion.
"""
import numpy as np
import yaml
from pathlib import Path
from ._utils import E_OVER_K

try:
    from numba import njit
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    def njit(func):
        return func

__all__ = ["ML_CVD_FDM"]


@njit
def _thomas(a, b, c_mat, rhs, cp, dp, out, n):
    cp[0] = c_mat[0] / b[0]
    dp[0] = rhs[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]
        cp[i] = c_mat[i] / m
        dp[i] = (rhs[i] - a[i] * dp[i - 1]) / m
    out[n - 1] = dp[n - 1]
    for i in range(n - 2, -1, -1):
        out[i] = dp[i] - cp[i] * out[i + 1]


@njit
def _interp_temp(t, t_arr, T_arr):
    n = len(t_arr)
    if t <= t_arr[0]:
        return T_arr[0]
    if t >= t_arr[-1]:
        return T_arr[-1]
    idx = 0
    for i in range(n - 1):
        if t_arr[i] <= t < t_arr[i + 1]:
            idx = i
            break
    f = (t - t_arr[idx]) / (t_arr[idx + 1] - t_arr[idx])
    return T_arr[idx] + f * (T_arr[idx + 1] - T_arr[idx])


@njit
def _compute_D(temp, d_coff, d_exp_c, d_temp_ref, dcal_type):
    if dcal_type == 0:
        return d_coff * 1e8 * np.exp(-d_exp_c * E_OVER_K / temp)
    else:
        return d_coff * np.exp(d_exp_c * (1.0 / temp - 1.0 / d_temp_ref))


@njit
def _evolve(
    c1_arr, c2_arr, x_arr, n_total,
    c1_init, c2_init,
    layer, pos, dep_t, n_layers,
    step_time, step_temperature,
    d1_coff, d1_exp_c, d1_temp_ref,
    d2_coff, d2_exp_c, d2_temp_ref,
    dcal_type, dx,
):
    n_active = 0
    dt_hist = np.zeros(500000)
    n_steps = 0
    D_hist = np.zeros((n_layers, 2))

    a_buf = np.zeros(n_total)
    b_buf = np.zeros(n_total)
    c_buf = np.zeros(n_total)
    rhs_buf = np.zeros(n_total)
    cp_buf = np.zeros(n_total)
    dp_buf = np.zeros(n_total)

    for ilayer in range(n_layers):
        thick = layer[ilayer]
        t_start = dep_t[ilayer]
        t_end = dep_t[ilayer + 1] if ilayer + 1 < len(dep_t) else dep_t[-1]
        t_rem = t_end - t_start
        if t_rem <= 0:
            continue

        if ilayer == 0:
            n_new = max(2, int(thick / dx))
            for j in range(n_new):
                x_arr[j] = pos[0] + j * dx
                c1_arr[j] = c1_init[0]
                c2_arr[j] = c2_init[0]
            n_active = n_new
            rate = 0.0
        else:
            rate = thick / t_rem

        t_elapsed = 0.0
        accum = 0.0

        while t_elapsed < t_rem - 1e-14:
            t_now = t_start + t_elapsed
            temp = _interp_temp(t_now, step_time, step_temperature) + 273.15

            D1 = _compute_D(temp, d1_coff, d1_exp_c, d1_temp_ref, dcal_type)
            D2 = _compute_D(temp, d2_coff, d2_exp_c, d2_temp_ref, dcal_type)
            D_max = max(abs(D1), abs(D2), 1e-30)
            dt_diff = 0.5 * dx * dx / D_max

            if ilayer == 0:
                dt = dt_diff
            else:
                dt_g = dx / rate * 1.5 if rate > 0 else 1.0
                dt = min(dt_g, dt_diff, 0.1)
            dt = max(dt, 1e-10)
            rem = t_rem - t_elapsed
            if dt > rem:
                dt = rem

            n_before = n_active

            if ilayer > 0 and rate > 0:
                accum += rate * dt
                n_new = int(accum / dx)
                if n_new > 0:
                    accum -= n_new * dx
                    end = min(n_active + n_new, n_total)
                    base = x_arr[n_active - 1]
                    for j in range(n_active, end):
                        x_arr[j] = base + (j - n_active + 1) * dx
                        c1_arr[j] = c1_init[ilayer]
                        c2_arr[j] = c2_init[ilayer]
                    n_active = end

            if n_active < 2:
                t_elapsed += dt
                continue

            r_val = D_max * dt / (2.0 * dx * dx)
            if r_val > 5.0:
                dt *= 0.3
                n_active = n_before
                continue

            n = n_active

            for i in range(n):
                a_buf[i] = -r_val
                b_buf[i] = 1.0 + 2.0 * r_val
                c_buf[i] = -r_val
            a_buf[0] = 0.0
            b_buf[0] = 1.0 + 2.0 * r_val
            c_buf[0] = -2.0 * r_val
            a_buf[n - 1] = -2.0 * r_val
            b_buf[n - 1] = 1.0 + 2.0 * r_val
            c_buf[n - 1] = 0.0

            for i in range(1, n - 1):
                rhs_buf[i] = (1.0 - 2.0 * r_val) * c1_arr[i] + r_val * (c1_arr[i - 1] + c1_arr[i + 1])
            rhs_buf[0] = c1_arr[0] * (1.0 - 2.0 * r_val) + 2.0 * r_val * c1_arr[1]
            rhs_buf[n - 1] = c1_arr[n - 1] * (1.0 - 2.0 * r_val) + 2.0 * r_val * c1_arr[n - 2]
            _thomas(a_buf, b_buf, c_buf, rhs_buf, cp_buf, dp_buf, c1_arr, n)

            for i in range(1, n - 1):
                rhs_buf[i] = (1.0 - 2.0 * r_val) * c2_arr[i] + r_val * (c2_arr[i - 1] + c2_arr[i + 1])
            rhs_buf[0] = c2_arr[0] * (1.0 - 2.0 * r_val) + 2.0 * r_val * c2_arr[1]
            rhs_buf[n - 1] = c2_arr[n - 1] * (1.0 - 2.0 * r_val) + 2.0 * r_val * c2_arr[n - 2]
            _thomas(a_buf, b_buf, c_buf, rhs_buf, cp_buf, dp_buf, c2_arr, n)

            if n_steps < len(dt_hist):
                dt_hist[n_steps] = dt
            n_steps += 1
            t_elapsed += dt

            if r_val < 0.01 and dt < 1.0:
                dt = min(dt * 1.2, 1.0)

        d1_s, d2_s = 0.0, 0.0
        tr = t_start
        dti = 1e-3
        while tr < t_end:
            tr_temp = _interp_temp(tr, step_time, step_temperature) + 273.15
            D1r = _compute_D(tr_temp, d1_coff, d1_exp_c, d1_temp_ref, dcal_type)
            D2r = _compute_D(tr_temp, d2_coff, d2_exp_c, d2_temp_ref, dcal_type)
            st = dti if tr + dti <= t_end else t_end - tr
            d1_s += D1r * st
            d2_s += D2r * st
            tr += st
        D_hist[ilayer, 0] = d1_s
        D_hist[ilayer, 1] = d2_s

    return c1_arr, c2_arr, x_arr, n_active, dt_hist, n_steps, D_hist


class ML_CVD_FDM:
    def __init__(self):
        self.c1_res = np.array([])
        self.c2_res = np.array([])
        self.x_position = np.array([])
        self.D = []
        self.dt_used = []

    def __call__(self, filepath: str, dcal_type: int = 0, dx: float = 1e-3) -> None:
        file_y = Path(filepath)
        with open(file_y, "r", encoding="utf-8") as fr:
            para_data = yaml.load(fr, Loader=yaml.FullLoader)

        layer = np.array(para_data["layer"], dtype=np.float64)
        dep_t = np.array(para_data["dep_t"], dtype=np.float64)
        c1_init = np.array(para_data["c1"], dtype=np.float64)
        c2_init = np.array(para_data["c2"], dtype=np.float64)

        d1_coff = para_data["d1_coff"]
        d1_temp_ref = para_data["d1_temp_ref"] + 273.15
        d1_exp_c = para_data["d1_exp_c"]
        d2_coff = para_data["d2_coff"]
        d2_temp_ref = para_data["d2_temp_ref"] + 273.15
        d2_exp_c = para_data["d2_exp_c"]

        step_temperature = np.array(para_data["step_temperature"], dtype=np.float64)
        step_time = np.array(para_data["step_time"], dtype=np.float64)

        pos = np.array([layer[:i].sum() for i in range(len(layer) + 1)])
        pos = pos - pos[1]

        n_layers = len(layer)
        total_thickness = pos[-1] - pos[0]
        n_total = int(np.ceil(total_thickness / dx)) + 2

        c1_arr = np.zeros(n_total)
        c2_arr = np.zeros(n_total)
        x_arr = np.zeros(n_total)

        c1_res, c2_res, x_res, n_active, dt_hist, n_steps, D_hist = _evolve(
            c1_arr, c2_arr, x_arr, n_total,
            c1_init, c2_init,
            layer, pos, dep_t, n_layers,
            step_time, step_temperature,
            d1_coff, d1_exp_c, d1_temp_ref,
            d2_coff, d2_exp_c, d2_temp_ref,
            dcal_type, dx,
        )

        self.c1_res = c1_res[:n_active].copy()
        self.c2_res = c2_res[:n_active].copy()
        self.x_position = x_res[:n_active].copy()
        self.D = [[D_hist[i, 0], D_hist[i, 1]] for i in range(n_layers)]
        self.dt_used = list(dt_hist[:n_steps])
