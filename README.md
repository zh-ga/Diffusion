# diffusion_sicm

**Multilayer Chemical Vapor Deposition (CVD) Diffusion Simulation**

[![CI](https://github.com/zh-ga/Diffusion/actions/workflows/build.yml/badge.svg)](https://github.com/zh-ga/Diffusion/actions/workflows/build.yml)

A Python library for simulating dopant diffusion in multilayer semiconductor structures during CVD processes. The default `ML_CVD_Model` solver uses a Crank-Nicolson finite difference method with Numba JIT acceleration. The original analytical erf solver is available as `ML_CVD_Model_Erf` for reference.

---

## Features

- **Default FDM solver** (`ML_CVD_Model`) — Crank-Nicolson with progressive grid growth
- **Backward compatible** — Same class name, same `__call__` signature, same output attributes
- **Sequential layer deposition** — Grid grows progressively as each layer is deposited
- **Adaptive time stepping** — Time step automatically adjusts to diffusion rate and growth rate
- **Numba JIT acceleration** — Inner loops compiled for near-native performance
- **Two-species diffusion** — Independent diffusion coefficients, activation energies, and concentration profiles
- **Arbitrary thermal history** — Temperature profile defined by arbitrary time-temperature curve
- **YAML configuration** — All process parameters in a single YAML file

---

## Installation

### Quick install from GitHub Release

```bash
pip install https://github.com/zh-ga/Diffusion/releases/download/v0.5.1/diffusion_sicm-0.5.1-py3-none-any.whl
```

### Install from source

```bash
git clone https://github.com/zh-ga/Diffusion.git
cd Diffusion
pip install -e .
```

---

## Quick Start

### 1. Prepare a YAML configuration file

Create `para_2layer.yaml`:

```yaml
layer: [5, 3]                    # layer thicknesses (um)
c1: [1.0e19, 1.0e15]            # species 1 conc. (cm^-3)
c2: [5.0e19, 1.0e15]            # species 2 conc. (cm^-3)
dep_t: [0, 100]                  # deposition time (s)

d1_coff: 0.76                    # D0 for species 1 (cm^2/s)
d1_temp_ref: 1100                # reference temperature (C)
d1_exp_c: 3.46                   # activation energy (eV)

d2_coff: 3.85                    # D0 for species 2 (cm^2/s)
d2_temp_ref: 1100
d2_exp_c: 3.66

step_temperature: [1100, 1100]    # temperature profile (C)
step_time: [0, 100]               # corresponding time points (s)
```

*Diffusion parameters for Boron (species 1) and Phosphorus (species 2) in silicon are based on published literature values.*

### 2. Run solver

```python
import matplotlib.pyplot as plt
from diffusion import ML_CVD_Model

model = ML_CVD_Model()
model("para_2layer.yaml", dcal_type=0)

plt.plot(model.x_position, model.c2_res, "r-")
plt.yscale("log")
plt.xlabel("Depth (um)")
plt.ylabel("Concentration (cm^-3)")
plt.show()
```

### 3. Legacy erf solver (for reference)

```python
from diffusion import ML_CVD_Model_Erf
model = ML_CVD_Model_Erf()
model("para_2layer.yaml", dcal_type=0)
```

---

## API Reference

### `ML_CVD_Model`

**Default solver** (FDM-based, replaces old analytical erf solver).

```python
model = ML_CVD_Model()
model(filepath, dcal_type=0, dx=1e-3)
```

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `filepath` | `str` | — | Path to YAML configuration file |
| `dcal_type` | `int` | `0` | Diffusion calculation mode: `0` = Arrhenius, `1` = exponential |
| `dx` | `float` | `1e-3` | Grid spacing in um (default 1 nm) |

**Attributes**

| Attribute | Type | Description |
|-----------|------|-------------|
| `c1_res` | `ndarray` | Species 1 concentration profile (cm^-3) |
| `c2_res` | `ndarray` | Species 2 concentration profile (cm^-3) |
| `x_position` | `ndarray` | Spatial coordinate (um) |
| `D` | `list` | Per-layer integrated diffusivities `[[D1_0, D2_0], ...]` |

---

### `ML_CVD_FDM`

Improved finite difference solver using Crank-Nicolson scheme with Numba JIT.

```python
model = ML_CVD_FDM()
model(filepath, dcal_type=0, dx=1e-3)
```

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `filepath` | `str` | — | Path to YAML configuration file |
| `dcal_type` | `int` | `0` | Diffusion calculation mode |
| `dx` | `float` | `1e-3` | Grid spacing in um (e.g. `1e-3` = 1 nm) |

**Attributes**

Same as `ML_CVD_Model`, plus:

| Attribute | Type | Description |
|-----------|------|-------------|
| `dt_used` | `list` | History of adaptive time steps used |

---

## Configuration File Format

```yaml
layer: [float, ...]           # thickness of each layer [um]
dep_t: [float, ...]           # deposition time points [s]
c1: [float, ...]              # species 1 concentration boundaries [cm^-3]
c2: [float, ...]              # species 2 concentration boundaries [cm^-3]

d1_coff: float                # D0 for species 1 [cm^2/s]
d1_temp_ref: float            # reference temperature [C]
d1_exp_c: float               # activation energy [eV] or exponent

d2_coff: float                # D0 for species 2
d2_temp_ref: float
d2_exp_c: float

step_temperature: [float, ...]  # temperature at each time point [C]
step_time: [float, ...]         # corresponding time points [s]
```

### Layer Interpretation

For `layer = [t0, t1, t2, t3]` (N entries):

| Internal Index | Physical Layer | Depth Range | Thickness |
|:---:|---|---|:---:|
| 0 | Substrate (base wafer) | `[-t0, 0]` | t0 |
| 1 | Layer 1 | `[0, t1]` | t1 |
| 2 | Layer 2 | `[t1, t1+t2]` | t2 |
| 3 | Layer 3 | `[t1+t2, t1+t2+t3]` | t3 |

Concentration arrays `c1`, `c2` have N entries, one per layer.

### Diffusion Calculation Modes

**Mode 0 — Arrhenius form:**

`D(T) = D0 x 10^8 x exp(-Ea x e / (k x T))`

**Mode 1 — Exponential form (relative to reference):**

`D(T) = D0 x exp(Ea x (1/T - 1/Tref))`

---

## Solver Comparison

| Feature | `ML_CVD_Model` (default, FDM) | `ML_CVD_Model_Erf` (legacy) |
|---------|:---:|:---:|
| Method | Crank-Nicolson FDM | Analytical erf solution |
| Status | **Recommended** | Legacy reference only |
| Layer coupling | Full (entire grid evolves together) | None (fixed boundary concentrations) |
| Deposition model | Progressive grid growth | Instantaneous full-layer |
| Grid adaptive | Yes (adaptive dt + r-value control) | No (fixed dx) |
| JIT compiled | Yes (Numba) | No |
| Physical accuracy | High | Moderate |

### Key Difference

The erf solver treats each layer independently with fixed boundary concentrations at interfaces. In contrast, the FDM solver solves the full PDE across the entire structure, allowing concentration at every point (including interfaces) to evolve continuously.

---

## Mathematical Background

### Governing Equation

Fick's second law for one-dimensional diffusion:

`dc/dt = D(T(t)) x d^2c/dx^2`

### Analytical Solution (erf)

For a diffusion couple with constant boundary concentrations cL and cR:

`c(x,t) = (cL + cR)/2 - (cL - cR)/2 x erf((x - x0) / (2 x sqrt(integral(D(t) dt))))`

### Numerical Solution (FDM)

Crank-Nicolson discretization:

```
(c_i^(n+1) - c_i^n) / dt = D/2 x (
  (c_(i-1)^(n+1) - 2c_i^(n+1) + c_(i+1)^(n+1)) / dx^2
  + (c_(i-1)^n - 2c_i^n + c_(i+1)^n) / dx^2
)
```

with Neumann (zero-flux) boundary conditions at both ends.

---

## Version Compatibility

| diffusion_sicm | Python  | Numba        | NumPy       | SciPy       |
|----------------|---------|--------------|-------------|-------------|
| 0.5.1          | >= 3.10 | >= 0.55      | >= 1.21     | >= 1.7      |

---

## Project Structure

```
diffusion_sicm/
├── LICENSE
├── README.md
├── CONTRIBUTING.md          # Contribution guide
├── pyproject.toml           # Build config + deps
├── .gitignore
├── .github/
│   └── workflows/
│       └── build.yml        # CI: test on push/PR, release on tag
├── examples/
│   ├── para_2layer.yaml     # Example: two-layer B/P diffusion in Si
│   ├── para_3layer.yaml     # Example: three-layer config
│   ├── run_erf.py
│   └── run_fdm.py
├── tests/
│   └── test_diffusion.py    # 15+ tests covering API, physics, edge cases
└── src/
    └── diffusion/
        ├── __init__.py       # Exports ML_CVD_Model, ML_CVD_FDM, ...
        ├── _version.py       # v0.5.1
        ├── _core.py          # ML_CVD_Model (FDM wrapper, backward compat)
        ├── _erf.py           # ML_CVD_Model_Erf (legacy erf solver)
        ├── _utils.py         # Shared utilities (interpolation, D calculation)
        └── fdm.py            # ML_CVD_FDM (Numba JIT Crank-Nicolson)
```

---

## License

Private — All rights reserved.

## Author

zh-ga
