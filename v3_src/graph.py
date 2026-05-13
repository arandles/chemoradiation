#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

dt_per_day = 24 * 60 * 2  # 2880 timesteps/day

base = Path(".")

sub_file = base / "sub" / "summary_fractional_change.csv"
opt_file = base / "opt" / "summary_fractional_change.csv"

sub = np.genfromtxt(sub_file, delimiter=",", names=True)
opt = np.genfromtxt(opt_file, delimiter=",", names=True)

sub_days = sub["timestep"] / dt_per_day
opt_days = opt["timestep"] / dt_per_day

plt.figure()

plt.plot(
    sub_days,
    sub["mean_fractional_change"],
    linewidth=2,
    label="Suboptimal mean"
)

plt.plot(
    opt_days,
    opt["mean_fractional_change"],
    linewidth=2,
    label="Optimal mean"
)

plt.fill_between(
    sub_days,
    sub["p25_fractional_change"],
    sub["p75_fractional_change"],
    alpha=0.25,
    label="Suboptimal IQR"
)

plt.fill_between(
    opt_days,
    opt["p25_fractional_change"],
    opt["p75_fractional_change"],
    alpha=0.25,
    label="Optimal IQR"
)

plt.xlabel("Days since therapy start")
plt.ylabel("Fractional volume change relative to therapy start")
plt.title("Tumor Volume Change")

plt.grid(True)

plt.xlim([0, max(np.max(sub_days), np.max(opt_days))])
plt.ylim([-1, 5])

plt.legend()

out_fig = base / "figure_3a_reproduction.png"

plt.savefig(out_fig, dpi=300, bbox_inches="tight")

plt.show()

print(f"Wrote graph to: {out_fig.resolve()}")
