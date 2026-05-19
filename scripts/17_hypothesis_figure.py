#!/usr/bin/env python3
"""
17_hypothesis_figure.py
Conceptual summary figure for the Discussion:
EtS skin microbiome equilibrium, dysbiosis trigger, and dual-pathway
hypothesis leading to equine summer eczema (sweet itch).

Run from project root:  python scripts/17_hypothesis_figure.py
Output: results/figures/discussion_eczema_hypothesis.pdf / .png
"""

import matplotlib
matplotlib.use("Agg")   # headless — no display required
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import os

os.makedirs("results/figures", exist_ok=True)

# ── Canvas ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(13, 16))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")
fig.patch.set_facecolor("white")

# ── Palette ───────────────────────────────────────────────────────────────────
C = dict(
    norm_bg   = "#c8e6c9",  # equilibrium background
    norm_box  = "#e8f5e9",  # sub-box inside equilibrium
    norm_edge = "#388e3c",
    trig_bg   = "#fff8e1",  # trigger
    trig_edge = "#f57f17",
    p1_bg     = "#fff3e0",  # pathway 1
    p2_bg     = "#fce4ec",  # pathway 2
    path_edge = "#c62828",
    bar_bg    = "#ef9a9a",  # barrier failure
    bar_edge  = "#b71c1c",
    ecz_bg    = "#b71c1c",  # eczema
    ecz_edge  = "#7f0000",
    arrow     = "#37474f",
    ref_bg    = "#eceff1",
    ref_edge  = "#90a4ae",
    text_dk   = "#212121",
    text_green= "#1b5e20",
    text_red  = "#b71c1c",
)

# ── Helpers ───────────────────────────────────────────────────────────────────
def box(x, y, w, h, fc, ec, lw=1.8, pad=0.025, zorder=2):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad={pad}",
        facecolor=fc, edgecolor=ec, linewidth=lw,
        transform=ax.transAxes, zorder=zorder, clip_on=False,
    ))

def arr(x0, y0, x1, y1, rad=0.0, lw=2.0, ms=18):
    ax.annotate(
        "", xy=(x1, y1), xytext=(x0, y0),
        xycoords="axes fraction", textcoords="axes fraction",
        arrowprops=dict(
            arrowstyle="->", color=C["arrow"],
            lw=lw, mutation_scale=ms,
            connectionstyle=f"arc3,rad={rad}",
        ), zorder=6,
    )

def txt(x, y, s, sz=8.5, w="normal", col="#212121",
        ha="center", va="center", style="normal", zorder=7):
    ax.text(x, y, s, transform=ax.transAxes,
            fontsize=sz, fontweight=w, fontstyle=style,
            color=col, ha=ha, va=va, zorder=zorder,
            linespacing=1.4)

# ═════════════════════════════════════════════════════════════════════════════
# TITLE
# ═════════════════════════════════════════════════════════════════════════════
txt(0.50, 0.977,
    "Elevated skin site microbiome and dual-pathway hypothesis\n"
    "for equine summer eczema",
    sz=12, w="bold", col="#1a237e")

# ═════════════════════════════════════════════════════════════════════════════
# TIER 1 — NORMAL EtS EQUILIBRIUM
# ═════════════════════════════════════════════════════════════════════════════
box(0.02, 0.848, 0.96, 0.108, C["norm_bg"], C["norm_edge"], lw=2.2)
txt(0.50, 0.950,
    "ELEVATED SKIN SITES (EtS) — normal microbial equilibrium",
    sz=10.5, w="bold", col=C["text_green"])

# Deinococcus sub-box
box(0.04, 0.854, 0.435, 0.086, C["norm_box"], C["norm_edge"])
txt(0.262, 0.927, "Deinococcaceae  (Deinococcus spp.)",
    sz=9, w="bold", col=C["text_green"])
txt(0.262, 0.908, "UV radiation & reactive oxygen species resistance",
    sz=8.2, col="#2e7d32")
txt(0.262, 0.891, "Carotenoid-mediated photoprotection  [1]",
    sz=8, style="italic", col="#388e3c")

# Micrococcus sub-box
box(0.525, 0.854, 0.455, 0.086, C["norm_box"], C["norm_edge"])
txt(0.752, 0.927, "Micrococcaceae  (Micrococcus spp.)",
    sz=9, w="bold", col=C["text_green"])
txt(0.752, 0.908, "Urea catabolism — balanced;  skin pH ≈ 5",
    sz=8.2, col="#2e7d32")
txt(0.752, 0.891, "Acid mantle maintenance  [2]",
    sz=8, style="italic", col="#388e3c")

# Arrow: equilibrium → trigger
arr(0.50, 0.848, 0.50, 0.791, lw=2.2)
txt(0.605, 0.820, "community\ndisruption",
    sz=7.5, style="italic", col=C["arrow"], ha="left")

# ═════════════════════════════════════════════════════════════════════════════
# DYSBIOSIS TRIGGER
# ═════════════════════════════════════════════════════════════════════════════
box(0.12, 0.718, 0.76, 0.073, C["trig_bg"], C["trig_edge"], lw=2.2)
txt(0.50, 0.772, "DYSBIOSIS TRIGGER",
    sz=11, w="bold", col="#e65100")
txt(0.50, 0.750,
    "Culicoides spp. bite  →  IgE-mediated hypersensitivity  →  pro-inflammatory cytokines  [3]",
    sz=8.8, col="#bf360c")

# Fork arrows: trigger → pathway 1 / pathway 2
arr(0.33, 0.718, 0.215, 0.658, rad=-0.15, lw=2.0)
arr(0.67, 0.718, 0.785, 0.658, rad=0.15,  lw=2.0)

# ═════════════════════════════════════════════════════════════════════════════
# PATHWAY 1 — Loss of Deinococcus
# ═════════════════════════════════════════════════════════════════════════════
box(0.01, 0.432, 0.455, 0.226, C["p1_bg"], C["path_edge"], lw=1.8)

txt(0.235, 0.650, "PATHWAY 1",
    sz=10, w="bold", col="#bf360c")
txt(0.235, 0.631, "Loss of Deinococcus",
    sz=9, w="bold", col=C["path_edge"])
txt(0.235, 0.612, "↓  UV photoprotection",
    sz=8.5, col=C["text_dk"])
txt(0.235, 0.594, "↓  Reactive oxygen species scavenging",
    sz=8.5, col=C["text_dk"])
txt(0.235, 0.576, "↓  Antioxidant defence capacity",
    sz=8.5, col=C["text_dk"])
txt(0.235, 0.553, "→  Oxidative damage to keratinocytes",
    sz=8.5, w="bold", col=C["text_red"])
txt(0.235, 0.535, "→  Impaired epidermal repair",
    sz=8.5, col=C["text_red"])
txt(0.235, 0.517, "→  Amplified local inflammatory response",
    sz=8.5, col=C["text_red"])
txt(0.235, 0.494, "[1]", sz=7.5, style="italic", col="#78909c")

# ═════════════════════════════════════════════════════════════════════════════
# PATHWAY 2 — Urea dysbiosis
# ═════════════════════════════════════════════════════════════════════════════
box(0.535, 0.432, 0.455, 0.226, C["p2_bg"], C["path_edge"], lw=1.8)

txt(0.762, 0.650, "PATHWAY 2",
    sz=10, w="bold", col="#bf360c")
txt(0.762, 0.631, "Micrococcus overgrowth / dysbiosis",
    sz=9, w="bold", col=C["path_edge"])
txt(0.762, 0.612, "↑  Urea hydrolysis  →  NH₃ + CO₂",
    sz=8.5, col=C["text_dk"])
txt(0.762, 0.594, "↑  Ammonia accumulation in sweat film",
    sz=8.5, col=C["text_dk"])
txt(0.762, 0.576, "↑  Skin surface pH  (> 6.5)",
    sz=8.5, col=C["text_dk"])
txt(0.762, 0.553, "→  Acid mantle disruption  [2, 4]",
    sz=8.5, w="bold", col=C["text_red"])
txt(0.762, 0.535, "→  Antimicrobial peptide inactivation",
    sz=8.5, col=C["text_red"])
txt(0.762, 0.517, "→  Compromised epidermal barrier",
    sz=8.5, col=C["text_red"])
txt(0.762, 0.494, "[2, 4]", sz=7.5, style="italic", col="#78909c")

# "+" additive symbol between the two pathways
txt(0.50, 0.545, "+", sz=28, w="bold", col=C["path_edge"])

# Merge arrows: pathways → barrier failure
arr(0.235, 0.432, 0.375, 0.373, rad=0.15, lw=2.0)
arr(0.762, 0.432, 0.625, 0.373, rad=-0.15, lw=2.0)

# ═════════════════════════════════════════════════════════════════════════════
# SKIN BARRIER FAILURE
# ═════════════════════════════════════════════════════════════════════════════
box(0.12, 0.285, 0.76, 0.088, C["bar_bg"], C["bar_edge"], lw=2.2)
txt(0.50, 0.349, "SKIN BARRIER FAILURE",
    sz=10.5, w="bold", col="white")
txt(0.50, 0.330,
    "Chronic pruritus  ·  Mast cell / IgE activation  ·  Secondary microbial colonisation",
    sz=8.8, col="#ffebee")
txt(0.50, 0.311, "Self-perpetuating itch–scratch cycle",
    sz=8.2, style="italic", col="#ffebee")

# Arrow: barrier → eczema
arr(0.50, 0.285, 0.50, 0.233, lw=2.5, ms=20)

# ═════════════════════════════════════════════════════════════════════════════
# EQUINE SUMMER ECZEMA
# ═════════════════════════════════════════════════════════════════════════════
box(0.10, 0.148, 0.80, 0.085, C["ecz_bg"], C["ecz_edge"], lw=2.8, zorder=3)
txt(0.50, 0.212,
    "EQUINE SUMMER ECZEMA  (Sweet itch)  [3]",
    sz=12, w="bold", col="white", zorder=8)
txt(0.50, 0.191,
    "Recurrent seasonal pruritic dermatitis; predilection for mane, tail base, and dorsum",
    sz=8.8, col="#ffcdd2", zorder=8)
txt(0.50, 0.172,
    "Clinical signs: alopecia, excoriation, lichenification, secondary infection",
    sz=8.2, style="italic", col="#ffcdd2", zorder=8)

# ═════════════════════════════════════════════════════════════════════════════
# REFERENCE BLOCK
# ═════════════════════════════════════════════════════════════════════════════
refs = (
    "[1] Battista JR (1997) Annu Rev Microbiol 51:203–224 — Deinococcus UV/ROS resistance\n"
    "[2] Kloos WE & Musselwhite MS (1975) Appl Microbiol 31:381–385 — Micrococcus skin ecology\n"
    "     Fluhr JW & Darlenski R (2014) Curr Probl Dermatol 49:1–10 — urea, pH and skin barrier\n"
    "[3] Schaffartzik A et al. (2012) Vet Immunol Immunopathol 147:113–126 — Culicoides hypersensitivity & sweet itch\n"
    "[4] Visscher MO et al. (2000) Skin Pharmacol Appl Skin Physiol 13:140–149 — urea hydrolysis and skin pH"
)
ax.text(0.02, 0.133, refs,
        transform=ax.transAxes,
        fontsize=7, color="#546e7a",
        va="top", ha="left",
        fontfamily="monospace",
        bbox=dict(facecolor=C["ref_bg"], edgecolor=C["ref_edge"],
                  boxstyle="round,pad=0.45", linewidth=0.8),
        zorder=4, linespacing=1.5)

# ═════════════════════════════════════════════════════════════════════════════
# Save
# ═════════════════════════════════════════════════════════════════════════════
plt.tight_layout(pad=0.3)
for fmt in ("pdf", "png"):
    out = f"results/figures/discussion_eczema_hypothesis.{fmt}"
    plt.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved {out}")
plt.show()
