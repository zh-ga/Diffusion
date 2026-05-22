"""Generate SVG images for LaTeX formulas used in README."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

OUTPUT_DIR = "docs/equations"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def render_latex(formula, filename, fontsize=18, figsize=None):
    """Render a LaTeX formula to SVG file with publication-quality fonts."""
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Liberation Serif', 'Times New Roman', 'DejaVu Serif'],
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Liberation Serif',
        'mathtext.it': 'Liberation Serif:italic',
        'mathtext.bf': 'Liberation Serif:bold',
        'mathtext.cal': 'Liberation Serif:italic',
        'mathtext.sf': 'Liberation Serif',
        'mathtext.tt': 'Liberation Serif',
        'mathtext.fallback': 'cm',
    })
    fig = plt.figure(figsize=figsize or (max(len(formula)*0.06, 3), 0.6))
    fig.text(0.02, 0.5, f"${formula}$", fontsize=fontsize,
             ha='left', va='center', color='black')
    plt.axis('off')
    filepath = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(filepath, format='svg', bbox_inches='tight', pad_inches=0.1,
                transparent=True, dpi=300)
    plt.close()
    print(f"  Created: {filepath}")

# ============================================================
# 1. Diffusion Calculation Modes
# ============================================================
render_latex(
    r"D(T) = D_0 \times 10^{8} \exp\left(-\frac{E_a e}{\mathit{kT}}\right)",
    "arrhenius_formula.svg",
    fontsize=16, figsize=(5.5, 0.6)
)

render_latex(
    r"D(T) = D_0 \exp\left(E_a \left(\frac{1}{\mathit{T}} - \frac{1}{\mathit{T_{ref}}}\right)\right)",
    "exponential_formula.svg",
    fontsize=16, figsize=(5.5, 0.6)
)

# ============================================================
# 2. Governing Equation (Fick's Second Law)
# ============================================================
render_latex(
    r"\frac{\partial c}{\partial t} = \mathit{D(T(t))} \frac{\partial^2 c}{\partial x^2}",
    "fick_law.svg",
    fontsize=18, figsize=(4.5, 0.6)
)

# ============================================================
# 3. Analytical Solution (erf)
# ============================================================
render_latex(
    r"\mathit{c(x,t)} = \frac{c_L+c_R}{2} - \frac{c_L-c_R}{2} \operatorname{erf}\left(\frac{x-x_0}{2\sqrt{D_{int}}}\right)",
    "erf_solution.svg",
    fontsize=15, figsize=(7.0, 0.6)
)

# ============================================================
# 4. Numerical Solution (Crank-Nicolson)
# ============================================================
render_latex(
    r"\frac{c_i^{n+1} - c_i^n}{\Delta t} = \frac{\mathit{D}}{2} \left( \frac{c_{i-1}^{n+1} - 2c_i^{n+1} + c_{i+1}^{n+1}}{\Delta x^2} + \frac{c_{i-1}^{n} - 2c_i^{n} + c_{i+1}^{n}}{\Delta x^2} \right)",
    "crank_nicolson.svg",
    fontsize=13, figsize=(9.0, 0.6)
)

print(f"\nAll SVG files saved to {OUTPUT_DIR}/")
