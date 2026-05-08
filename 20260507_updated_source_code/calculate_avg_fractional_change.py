#!/usr/bin/env python3
import csv
from pathlib import Path
import numpy as np

def read_two_col_csv(csv_path: Path):
    """
    Reads a 2-col CSV: timestep,value
    Returns: dict {timestep:int -> value:float}
    Skips header / bad lines.
    """
    series = {}
    with csv_path.open("r", newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or len(row) < 2:
                continue
            try:
                t = int(str(row[0]).strip())
                v = float(str(row[1]).strip())
            except ValueError:
                continue
            series[t] = v
    if not series:
        raise ValueError(f"No valid data rows found in {csv_path}")
    return series

def main(
    base_dir=".",
    i_start=1,
    i_end=500,
    filename_template="output_all_{i}.csv",
    out_file="avg_fractional_change.csv",
):
    base = Path(base_dir)

    all_series = []
    used_i = []

    for i in range(i_start, i_end + 1):
        csv_path = base / filename_template.format(i=i)
        if not csv_path.exists():
            print(f"[WARN] Missing: {csv_path}")
            continue

        try:
            raw = read_two_col_csv(csv_path)
        except Exception as e:
            print(f"[WARN] Failed reading {csv_path}: {e}")
            continue

        # Normalize each replicate to fractional change from its own start:
        # frac_change(t) = v(t)/v0 - 1
        t0 = min(raw.keys())
        v0 = raw[t0]

        if v0 == 0:
            print(f"[WARN] Skipping {csv_path}: starting value is 0")
            continue

        frac = {t: (v / v0) - 1.0 for t, v in raw.items()}

        all_series.append(frac)
        used_i.append(i)

    if not all_series:
        raise SystemExit("No valid input files found. Check filenames and directory.")

    # Collect all values by timestep
    vals_by_t = {}
    for s in all_series:
        for t, v in s.items():
            vals_by_t.setdefault(t, []).append(v)

    timesteps = sorted(vals_by_t.keys())

    out_path = base / out_file
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "timestep",
            "avg_fractional_change",
            "lower_bound_p25",
            "upper_bound_p75",
            "n_files"
        ])

        for t in timesteps:
            vals = np.array(vals_by_t[t], dtype=float)
            avg = np.mean(vals)
            lower = np.percentile(vals, 25)
            upper = np.percentile(vals, 75)
            writer.writerow([t, avg, lower, upper, len(vals)])

    print(f"Used files i: {used_i}")
    print(f"Wrote: {out_path} (rows: {len(timesteps)})")

if __name__ == "__main__":
    main()

