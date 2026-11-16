"""Evidence plots: |B|, speed, density and temperature around each detection.

    python scripts/plot_candidates.py wind data/cache results/benchmark_..._detections.csv out.png [--which unmatched-passed]

One column per event, 30 minutes either side, with the detection time marked
and the upstream and downstream averaging windows shaded.
"""

from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from helios.data.cache import load_cached

BG, FG, MUTED, GRID = "#0b0e13", "#e9e6df", "#8b93a1", "#1e242e"
COLORS = {"b": "#e8b45a", "v": "#7fb2e5", "n": "#6fd0c0", "t": "#ef8a6b"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("spacecraft")
    ap.add_argument("cache")
    ap.add_argument("detections")
    ap.add_argument("out")
    ap.add_argument("--which", default="unmatched-passed", choices=["unmatched-passed", "unmatched-failed", "matched"])
    ap.add_argument("--max", type=int, default=14)
    args = ap.parse_args()

    data = load_cached(args.spacecraft, args.cache)
    t = pd.read_csv(args.detections, parse_dates=["time"])
    sel = {"unmatched-passed": (~t.matched) & t.audit_passed, "unmatched-failed": (~t.matched) & ~t.audit_passed,
           "matched": t.matched}[args.which]
    ev = t[sel].head(args.max).reset_index(drop=True)
    cols = 7
    rows = int(np.ceil(len(ev) / cols))
    fig, axes = plt.subplots(4 * rows, cols, figsize=(3.0 * cols, 5.2 * rows), squeeze=False, facecolor=BG)
    half, minute = pd.Timedelta(minutes=30), pd.Timedelta(minutes=1)
    for k, r in ev.iterrows():
        row0, c = 4 * (k // cols), k % cols
        m = data.mag.loc[r.time - half: r.time + half]
        p = data.plasma.loc[r.time - half: r.time + half]
        x = lambda idx: (idx - r.time).total_seconds() / 60  # noqa: E731
        series = [("b", m.index, m["bmag"], "|B| nT"), ("v", p.index, np.linalg.norm(p[["vx", "vy", "vz"]].to_numpy(), axis=1), "V km/s"),
                  ("n", p.index, p["n"], "N cm⁻³"), ("t", p.index, p["tp"] / 1e4, "T 10⁴ K")]
        up, down = ((-9, -1), (2, 10)) if r.type == "FF" else ((1, 9), (-10, -2))
        for i, (key, idx, y, label) in enumerate(series):
            ax = axes[row0 + i][c]
            ax.set_facecolor(BG)
            ax.axvspan(*up, color="#ffffff", alpha=0.05, lw=0)
            ax.axvspan(*down, color="#ffffff", alpha=0.10, lw=0)
            ax.axvline(0, color=FG, lw=0.6, alpha=0.6)
            ax.plot(x(idx), y, color=COLORS[key], lw=1.0)
            ax.set_xlim(-30, 30)
            ax.tick_params(colors=MUTED, labelsize=7, length=2)
            ax.grid(color=GRID, lw=0.5)
            for s in ax.spines.values():
                s.set_color(GRID)
            if c == 0:
                ax.set_ylabel(label, color=MUTED, fontsize=8)
            if i == 0:
                ax.set_title(f"{r.time:%Y-%m-%d %H:%M}  {r.type}\nM_f {r.mach_fast:.2f}  sharp {r.sharpness:.2f}  B×{r.b_ratio:.2f}",
                             color=FG, fontsize=8, loc="left")
            if i < 3:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel("minutes", color=MUTED, fontsize=7)
    for k in range(len(ev), rows * cols):
        for i in range(4):
            axes[4 * (k // cols) + i][k % cols].set_visible(False)
    fig.tight_layout(h_pad=0.3, w_pad=0.6)
    fig.savefig(args.out, dpi=110, facecolor=BG)
    print(f"wrote {args.out} with {len(ev)} events")


if __name__ == "__main__":
    main()
