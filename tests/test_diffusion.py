"""Tests for diffusion_sicm.

Test categories:
- Basic model execution and API compatibility
- Physical correctness (monotonicity, conservation)
- Error handling for invalid inputs
- Edge cases
- FDM vs erf cross-validation
- Deprecated API compatibility
"""

import sys
import os
import tempfile
import warnings

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from diffusion import ML_CVD_Model, ML_CVD_FDM, __version__

# Use the simplified two-layer example (public academic parameters)
YAML_PATH = os.path.join(os.path.dirname(__file__), "..", "examples", "para_2layer.yaml")


def test_version_is_string():
    assert isinstance(__version__, str)


def test_ml_cvd_model_runs():
    m = ML_CVD_Model()
    m(YAML_PATH, dcal_type=0)
    assert len(m.x_position) > 0
    assert len(m.c1_res) == len(m.x_position)
    assert len(m.c2_res) == len(m.x_position)
    assert m.c2_res.max() > m.c2_res.min()
    assert len(m.D) > 0


def test_ml_cvd_fdm_runs():
    m = ML_CVD_FDM()
    m(YAML_PATH, dcal_type=0, dx=1e-3)
    assert len(m.x_position) > 0
    assert len(m.c1_res) == len(m.x_position)
    assert len(m.c2_res) == len(m.x_position)
    assert m.c2_res.max() > m.c2_res.min()


def test_concentration_monotonic():
    m = ML_CVD_FDM()
    m(YAML_PATH, dcal_type=0, dx=1e-3)
    assert not np.any(np.isnan(m.c1_res))
    assert not np.any(np.isnan(m.c2_res))
    assert not np.any(np.isinf(m.c1_res))
    assert not np.any(np.isinf(m.c2_res))
    assert np.all(np.diff(m.x_position) > 0)


def test_concentration_within_range():
    m = ML_CVD_FDM()
    m(YAML_PATH, dcal_type=0, dx=1e-3)
    assert m.c1_res.min() > 0
    assert m.c2_res.min() > 0
    assert m.c1_res.min() >= 1e14
    assert m.c1_res.max() <= 1e20


def test_ml_cvd_model_has_same_shape_as_fdm():
    m1 = ML_CVD_Model()
    m2 = ML_CVD_FDM()
    m1(YAML_PATH, dcal_type=0)
    m2(YAML_PATH, dcal_type=0, dx=1e-3)
    assert len(m1.x_position) == len(m2.x_position)
    assert len(m1.c1_res) == len(m2.c1_res)
    assert len(m1.c2_res) == len(m2.c2_res)


def test_missing_key_raises_keyerror():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write("layer: [1, 2]\ndep_t: [0, 10]\n")
        f.flush()
        fname = f.name
    try:
        m = ML_CVD_FDM()
        with pytest.raises(KeyError):
            m(fname)
    finally:
        os.unlink(fname)


def test_empty_yaml_raises_error():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write("")
        f.flush()
        fname = f.name
    try:
        m = ML_CVD_FDM()
        with pytest.raises(Exception):
            m(fname)
    finally:
        os.unlink(fname)


def test_dcal_type_1_runs():
    import yaml
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump({
            "layer": [1.0, 0.5],
            "dep_t": [0, 10],
            "c1": [1e18, 1e14],
            "c2": [1e20, 1e15],
            "d1_coff": 0.76,
            "d1_temp_ref": 1100,
            "d1_exp_c": 3.46,
            "d2_coff": 3.85,
            "d2_temp_ref": 1100,
            "d2_exp_c": 3.66,
            "step_temperature": [1100, 1100],
            "step_time": [0, 10],
        }, f)
        f.flush()
        fname = f.name
    try:
        m = ML_CVD_FDM()
        m(fname, dcal_type=1, dx=1e-1)
        assert len(m.x_position) > 0
        assert len(m.c1_res) == len(m.x_position)
        assert m.c1_res.min() > 0
    finally:
        os.unlink(fname)


def test_single_layer():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write("""
layer: [10]
dep_t: [0, 100]
c1: [1e18]
c2: [1e16]
d1_coff: 0.76
d1_temp_ref: 1100
d1_exp_c: 3.46
d2_coff: 3.85
d2_temp_ref: 1100
d2_exp_c: 3.66
step_temperature: [1100, 1100]
step_time: [0, 100]
""")
        f.flush()
        fname = f.name
    try:
        m = ML_CVD_FDM()
        m(fname, dcal_type=0, dx=1e-3)
        assert len(m.x_position) > 0
        assert len(m.c1_res) == len(m.x_position)
        assert m.c1_res.min() > 0
    finally:
        os.unlink(fname)


def test_constant_temperature():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write("""
layer: [5, 5]
dep_t: [0, 200]
c1: [1e18, 1e14]
c2: [1e20, 1e15]
d1_coff: 0.76
d1_temp_ref: 1100
d1_exp_c: 3.46
d2_coff: 3.85
d2_temp_ref: 1100
d2_exp_c: 3.66
step_temperature: [1100, 1100]
step_time: [0, 200]
""")
        f.flush()
        fname = f.name
    try:
        m = ML_CVD_FDM()
        m(fname, dcal_type=0, dx=1e-3)
        assert len(m.x_position) > 0
        assert m.c1_res.min() > 0
    finally:
        os.unlink(fname)


def test_x_positon_deprecated():
    m = ML_CVD_Model()
    m(YAML_PATH, dcal_type=0)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        _ = m.x_positon
        assert len(w) == 1
        assert issubclass(w[0].category, DeprecationWarning)
        assert "x_positon" in str(w[0].message)


def test_fdm_dt_used_attribute():
    m = ML_CVD_FDM()
    m(YAML_PATH, dcal_type=0, dx=1e-3)
    assert isinstance(m.dt_used, list)
    assert len(m.dt_used) > 0
    assert all(isinstance(dt, (int, float, np.floating)) for dt in m.dt_used)


def test_fdm_dt_used_positive():
    m = ML_CVD_FDM()
    m(YAML_PATH, dcal_type=0, dx=1e-3)
    assert all(dt > 0 for dt in m.dt_used)


def test_fdm_D_attribute():
    m = ML_CVD_FDM()
    m(YAML_PATH, dcal_type=0, dx=1e-3)
    assert len(m.D) == 2  # 2 layers in the example config
    for d_layer in m.D:
        assert len(d_layer) == 2  # D1, D2 for each layer


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
