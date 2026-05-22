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
- **Numba JIT acceleration** — Inner loops compiled for near-native performance (~3s for 400s process at 1nm resolution)
- **Two-species diffusion** — Independent diffusion coefficients, activation energies, and concentration profiles
- **Arbitrary thermal history** — Temperature profile defined by arbitrary time-temperature curve
- **YAML configuration** — All process parameters in a single YAML file

---

## Installation

### Quick install from GitHub Release

```bash
pip install https://github.com/zh-ga/Diffusion/releases/download/v0.5.0/diffusion_sicm-0.5.0-py3-none-any.whl
```

### Install from source

```bash
git clone https://github.com/zh-ga/Diffusion.git
cd Diffusion
pip install -e .
# or with FDM extras:
pip install -e ".[fdm]"
```

---

## Quick Start

### 1. Prepare a YAML configuration file

```yaml
# para_3layer.yaml
layer: [4, 0.5, 6.2, 3.5]      # layer thicknesses (um)
c1: [1.01e10, 1.01e10, 1.01e10, 1.01e10]  # species 1 conc. (cm^-3)
c2: [3.48e19, 3.96e10, 2.88e16, 1.10e16]  # species 2 conc. (cm^-3)
dep_t: [0, 30, 250, 400]        # deposition times (s)

d1_coff: 2.38                   # diffusion pre-factor (species 1)
d1_temp_ref: 1120               # reference temperature (C)
d1_exp_c: 3.60                  # activation energy (eV, dcal_type=0) or exponent

d2_coff: 2.62                   # diffusion pre-factor (species 2)
d2_temp_ref: 1120
d2_exp_c: 3.62

step_temperature: [1110, 1110, 1110, 1110, 1110, 1110, 1110, 700]
step_time: [0, 20, 30, 230, 250, 350, 370, 400]
```

### 2. Run solver

```python
import matplotlib.pyplot as plt
from diffusion import ML_CVD_Model

model = ML_CVD_Model()
model("para_3layer.yaml", dcal_type=0)

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
model("para_3layer.yaml", dcal_type=0)
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
$$D(T) = D_0 \\times 10^{8} \\times \\exp\\left(-\\frac{E_a \\cdot e}{k \\cdot T}\\right)$$

**Mode 1 — Exponential form (relative to reference):**
$$D(T) = D_0 \\times \\exp\\left(E_a \\cdot \\left(\\frac{1}{T} - \\frac{1}{T_{ref}}\\right)\\right)$$

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
| Computation time (400s) | ~3s (warm) | ~0.08s |
| Physical accuracy | High | Moderate |

### Key Difference

The erf solver treats each layer independently with fixed boundary concentrations at interfaces. This means the concentration at a given interface **does not evolve** once set. In contrast, the FDM solver solves the full PDE across the entire structure, allowing concentration at every point (including interfaces) to evolve continuously in response to diffusion from neighboring layers.

---

## Benchmark Results

Tested on `para_3layer.yaml` (400s process, dx = 1 nm):

```
Metric                              Original erf       FDM v2
------------------------------------------------------------
Computation time (s)                      0.0835       2.9957
Grid points                               141999        10700
C2 min (Layer 1 & 2)                  2.9285e+14   7.5589e+15
C2 max (Layer 1 & 2)                  3.4791e+19   3.4791e+19
```

The large difference in C2 minimum (26x) is **physical** — the erf model's fixed low boundary concentration at the Layer 1/2 interface causes excessive depletion in Layer 1, while the FDM model correctly captures the continuous supply of dopant from Layer 2.

---

## Mathematical Background

### Governing Equation

Fick's second law for one-dimensional diffusion:

$$\\frac{\\partial c}{\\partial t} = D(T(t)) \\frac{\\partial^2 c}{\\partial x^2}$$

### Analytical Solution (erf)

For a diffusion couple with constant boundary concentrations $c_L$ and $c_R$:

$$c(x, t) = \\frac{c_L + c_R}{2} - \\frac{c_L - c_R}{2} \\cdot \\operatorname{erf}\\left(\\frac{x - x_0}{2\\sqrt{\\int D(t) dt}}\\right)$$

### Numerical Solution (FDM)

Crank-Nicolson discretization:

$$\\frac{c_i^{n+1} - c_i^n}{\\Delta t} = \\frac{D}{2} \\left( \\frac{c_{i-1}^{n+1} - 2c_i^{n+1} + c_{i+1}^{n+1}}{\\Delta x^2} + \\frac{c_{i-1}^{n} - 2c_i^{n} + c_{i+1}^{n}}{\\Delta x^2} \\right)$$

with Neumann (zero-flux) boundary conditions at both ends.

---

## Version Compatibility

| diffusion_sicm | Python  | Numba        | NumPy       | SciPy       |
|----------------|---------|--------------|-------------|-------------|
| 0.4.0          | >= 3.10 | >= 0.55      | >= 1.21     | >= 1.7      |

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
│   ├── para_3layer.yaml
│   ├── run_erf.py
│   └── run_fdm.py
├── tests/
│   └── test_diffusion.py    # 15+ tests covering API, physics, edge cases
└── src/
    └── diffusion/
        ├── __init__.py       # Exports ML_CVD_Model, ML_CVD_FDM, ...
        ├── _version.py       # v0.4.0
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
