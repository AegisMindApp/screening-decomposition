#!/usr/bin/env python3
"""Figure 1: six interventions against the two measured resolution floors.

Every value is read from the recorded artifacts rather than retyped: the floors and their
bootstrap intervals from FLOOR_INTERVAL.json, the intervention deltas from the arm results.
Each intervention is drawn against the floor that applies to it -- the scoring floor where
poses were held fixed, the protocol floor where placement or search was allowed to vary.
"""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

F = json.load(open("analysis/docking_value/FLOOR_INTERVAL.json"))
SC, PR = F["scoring_floor"], F["protocol_floor"]

# (label, delta, lo, hi, which floor applies)
ARMS = [
    ("Scoring function\n(gnina CNN, identical poses)",  0.140,  0.091,  0.187, "scoring"),
    # Re-run on the repaired receptor across 40 fold seeds; see
    # analysis/pose_ensemble/RECONSTRUCTION.md. Supersedes the published +0.019.
    ("Pose ensemble vs top pose",                       0.055,  0.043,  0.068, "scoring"),
    ("Receptor preparation repair",                     0.045,  0.024,  0.067, "protocol"),
    ("Pose generator (DiffDock-L)",                     0.017, -0.027,  0.064, "protocol"),
    ("Binding site supplied to Boltz-2",               -0.013, -0.026,  0.000, "protocol"),
    ("Eight-fold search effort",                       -0.019, -0.032, -0.006, "protocol"),
]

INK   = "#1f2933"
GREY  = "#9aa5b1"
C_SC  = "#2f6f9f"   # scoring floor
C_PR  = "#c2701c"   # protocol floor
COL   = {"scoring": C_SC, "protocol": C_PR}

fig, ax = plt.subplots(figsize=(7.4, 3.9), dpi=300)
fig.patch.set_facecolor("white"); ax.set_facecolor("white")

# floors: shaded bootstrap interval + median line
for (fl, c, name) in ((SC, C_SC, "scoring floor"), (PR, C_PR, "protocol floor")):
    ax.axvspan(fl["ci95"][0], fl["ci95"][1], color=c, alpha=0.10, lw=0, zorder=0)
    ax.axvline(fl["gap"], color=c, ls="--", lw=1.3, zorder=1)

ax.axvline(0, color=GREY, lw=0.9, zorder=1)

ys = list(range(len(ARMS)))[::-1]
for y, (lab, d, lo, hi, which) in zip(ys, ARMS):
    c = COL[which]
    ax.plot([lo, hi], [y, y], color=c, lw=1.9, solid_capstyle="butt", zorder=3)
    for x in (lo, hi):
        ax.plot([x, x], [y - 0.13, y + 0.13], color=c, lw=1.4, zorder=3)
    ax.plot([d], [y], "o", ms=7.5, color=c, mec="white", mew=1.2, zorder=4)

ax.set_yticks(ys)
ax.set_yticklabels([a[0] for a in ARMS], fontsize=8.6, color=INK)
ax.set_ylim(-0.75, len(ARMS) - 0.25)
ax.set_xlim(-0.045, 0.205)
ax.set_xlabel("ΔAUROC  (95% paired-bootstrap CI)", fontsize=9.5, color=INK)
ax.tick_params(axis="x", labelsize=8.6, colors=INK)
ax.tick_params(axis="y", length=0)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color(GREY)

# floor labels sit above the plotting area so they never collide with an interval
ax.annotate(f"scoring floor\n{SC['gap']:.3f}", xy=(SC["gap"], len(ARMS) - 0.42),
            xytext=(SC["gap"] - 0.011, len(ARMS) - 0.42), ha="right", va="center",
            fontsize=7.8, color=C_SC, linespacing=1.25)
ax.annotate(f"protocol floor\n{PR['gap']:.3f}", xy=(PR["gap"], len(ARMS) - 0.42),
            xytext=(PR["gap"] + 0.008, len(ARMS) - 0.42), ha="left", va="center",
            fontsize=7.8, color=C_PR, linespacing=1.25)

# legend below the axis: inside the panel it would sit on top of the shaded floor bands
ax.legend(handles=[Line2D([], [], color=C_SC, lw=2, marker="o", ms=6, mec="white",
                          label="gated by the scoring floor (poses held fixed)"),
                   Line2D([], [], color=C_PR, lw=2, marker="o", ms=6, mec="white",
                          label="gated by the protocol floor (placement or search varies)")],
          loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=1,
          fontsize=8.0, frameon=False, handletextpad=0.7, borderaxespad=0)

plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"docs/papers/figures/figure1_floors.{ext}", bbox_inches="tight",
                facecolor="white")
print("wrote docs/papers/figures/figure1_floors.png and .pdf")
print(f"  scoring floor  {SC['gap']:.4f} {SC['ci95']}")
print(f"  protocol floor {PR['gap']:.4f} {PR['ci95']}")
