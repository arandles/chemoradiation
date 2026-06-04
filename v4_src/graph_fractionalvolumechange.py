#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

dt_per_day = 24 * 60 * 2  # 2880 timesteps/day

base = Path(".")

sub_file  = base / "sub"  / "summary_fractional_change.csv"
opt_file  = base / "opt"  / "summary_fractional_change.csv"
soc_file  = base / "SOC"  / "summary_fractional_change.csv"
rand_file = base / "rand" / "summary_fractional_change.csv"

for f in [sub_file, opt_file, soc_file, rand_file]:
    if not f.exists():
        raise FileNotFoundError(f"Missing file: {f.resolve()}")

sub  = np.genfromtxt(sub_file,  delimiter=",", names=True)
opt  = np.genfromtxt(opt_file,  delimiter=",", names=True)
soc  = np.genfromtxt(soc_file,  delimiter=",", names=True)
rand = np.genfromtxt(rand_file, delimiter=",", names=True)

sub_days  = sub["timestep"]  / dt_per_day
opt_days  = opt["timestep"]  / dt_per_day
soc_days  = soc["timestep"]  / dt_per_day
rand_days = rand["timestep"] / dt_per_day

plt.figure(figsize=(8, 5))

plt.plot(
    sub_days,
    sub["mean_fractional_change"],
    linewidth=2,
    label="Suboptimal"
)

plt.plot(
    opt_days,
    opt["mean_fractional_change"],
    linewidth=2,
    label="Optimal"
)

plt.plot(
    soc_days,
    soc["mean_fractional_change"],
    linewidth=2,
    linestyle="--",
    marker="o",
    markevery=200,
    markersize=4,
    label="SOC"
)

plt.plot(
    rand_days,
    rand["mean_fractional_change"],
    linewidth=2,
    linestyle=":",
    marker="s",
    markevery=200,
    markersize=4,
    label="Random"
)

plt.axhline(0, linewidth=1, linestyle="--")

plt.xlabel("Time (days)")
plt.ylabel("Fractional change from initial value")
plt.title("Mean fractional change over time")
plt.ylim(-1, 3.5)
plt.legend()
plt.tight_layout()

plt.savefig("fractional_change_sub_opt_soc_rand_lines_only.png", dpi=300)
plt.savefig("fractional_change_sub_opt_soc_rand_lines_only.pdf")

plt.show()
