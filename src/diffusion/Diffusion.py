"""Backward-compatible module alias for existing user code.

This module preserves the original import path:
    from diffusion.Diffusion import ML_CVD_Model
so that existing scripts continue to work without changes.
"""

from ._core import ML_CVD_Model
from ._erf import ML_CVD_Model_Erf

__all__ = ["ML_CVD_Model", "ML_CVD_Model_Erf"]
