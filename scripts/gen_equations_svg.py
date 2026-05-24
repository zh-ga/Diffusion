"""Generate publication-quality formula cards for README.

Uses Computer Modern math fonts (standard LaTeX look) with
a clean, academic-style dark card design.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

OUTPUT_DIR = "docs/equations"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Color palette ──────────────────────────────────────────
BG_TOP   = "#1c1c2a"
BG_BOT   = "#14141e"
TEXT_CLR = "#e8e8e8"
ACCENT   = "#6fa8dc"     # muted blue, academic feel


def render_formula_card(
    formula: str,
    filename: str,
    fontsize: int = 18,
    figsize: tuple = (6.0, 0.8),
) -> None:
    """Render a LaTeX formula on an academic-style card as SVG."""
    w, h = figsize

    fig = plt.figure(figsize=(w, h), dpi=200)
    fig.patch.set_visible(False)

    # ── background gradient ────────────────────────────────
    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    ax_bg = fig.add_axes([0, 0, 1, 1], zorder=0)
    ax_bg.imshow(
        gradient,
        extent=[0, w, 0, h],
        aspect="auto",
        cmap=plt.matplotlib.colors.LinearSegmentedColormap.from_list(
            "bg", [BG_TOP, BG_BOT]
        ),
    )
    ax_bg.axis("off")

    # ── subtle accent line at top ──────────────────────────
    accent_line = mpatches.FancyBboxPatch(
        (0.02, h - 0.04),
        w - 0.04,
        0.02,
        boxstyle=mpatches.BoxStyle("Round", pad=0.005),
        linewidth=0,
        facecolor=ACCENT,
        alpha=0.6,
        transform=fig.dpi_scale_trans,
    )
    ax_bg.add_patch(accent_line)

    # ── formula text (Computer Modern for true LaTeX look) ─
    rc_params = {
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "DejaVu Serif", "Liberation Serif"],
        "mathtext.fontset": "cm",        # ← Computer Modern (LaTeX standard)
        "mathtext.default": "it",
        "text.color": TEXT_CLR,
        "axes.edgecolor": "none",
        "axes.facecolor": "none",
    }
    with plt.rc_context(rc_params):
        ax_text = fig.add_axes([0.05, 0.10, 0.90, 0.82], zorder=2)
        ax_text.text(
            0.5, 0.5,
            f"${formula}$",
            fontsize=fontsize,
            color=TEXT_CLR,
            ha="center",
            va="center",
            transform=ax_text.transAxes,
        )
        ax_text.axis("off")

    filepath = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(
        filepath,
        format="svg",
        bbox_inches="tight",
        pad_inches=0.10,
        facecolor=fig.get_facecolor(),
        transparent=False,
    )
    plt.close()
    print(f"  Created: {filepath}")


# ═══════════════════════════════════════════════════════════
#  Formulas – keep LaTeX simple for matplotlib's mathtext
# ═══════════════════════════════════════════════════════════

formulas = [
    ("arrhenius_formula.svg",
     r"D(T) = D_0 \times 10^{8} \exp\left(-\frac{E_a e}{kT}\right)",
     18, (5.5, 0.85)),

    ("exponential_formula.svg",
     r"D(T) = D_0 \exp\left(E_a \left(\frac{1}{T} - \frac{1}{T_{\mathrm{ref}}}\right)\right)",
     18, (5.5, 0.85)),

    ("fick_law.svg",
     r"\frac{\partial c}{\partial t} = D(T(t)) \frac{\partial^2 c}{\partial x^2}",
     20, (4.8, 0.85)),

    ("erf_solution.svg",
     r"c(x,t) = \frac{c_L + c_R}{2} - \frac{c_L - c_R}{2} \operatorname{erf}\left(\frac{x - x_0}{2\sqrt{D_{\mathrm{int}}}}\right)",
     16, (7.8, 0.85)),

    ("crank_nicolson.svg",
     r"\frac{c_i^{n+1} - c_i^n}{\Delta t} = \frac{D}{2} \left["
     r"\frac{c_{i-1}^{n+1} - 2c_i^{n+1} + c_{i+1}^{n+1}}{\Delta x^2}"
     r"+ \frac{c_{i-1}^{n} - 2c_i^{n} + c_{i+1}^{n}}{\Delta x^2} \right]",
     14, (10.5, 0.85)),
]

for fname, formula, fsize, fsize_fig in formulas:
    render_formula_card(formula, fname, fontsize=fsize, figsize=fsize_fig)

print(f"\n✅ All {len(formulas)} formula cards saved to {OUTPUT_DIR}/")
