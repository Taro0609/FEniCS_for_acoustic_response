"""
Extract target-point responses from FEniCSx/dolfinx VTU outputs.

Required input files
--------------------
target_coordinates.csv
    CSV file containing measurement points. Rows with the following
    point_name values are read by default:

        columella
        footplate
        prootic

pressure_p0_*.vtu
    Pressure outputs from frequency_sweep_solver.py.

displacement_p0_*.vtu
    Displacement outputs from frequency_sweep_solver.py.

Outputs
-------
target_response.csv
    Frequency-response values at each target point.

target_response_plot.png
    Summary plot of displacement response at target points.
"""

import argparse
import csv
from pathlib import Path

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt


def read_points_from_csv(csv_path, point_names):
    """Read named 3-D points from target_coordinates.csv."""
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Coordinate file not found: {csv_path}")

    wanted = set(point_names)
    points = {}

    with csv_path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get("point_name")
            if name in wanted:
                points[name] = np.array(
                    [float(row["x"]), float(row["y"]), float(row["z"])],
                    dtype=float,
                )

    missing = wanted - set(points)
    if missing:
        raise ValueError(f"Missing target points in {csv_path}: {sorted(missing)}")

    return points


def closest_point_index(mesh_points, target):
    distances = np.linalg.norm(mesh_points - target, axis=1)
    return int(np.argmin(distances))


def get_array(grid, preferred_names):
    for name in preferred_names:
        if name in grid.point_data:
            return grid.point_data[name]
    available = ", ".join(grid.point_data.keys())
    raise KeyError(f"None of {preferred_names} found. Available arrays: {available}")


def displacement_norm_at_file(filename, target):
    grid = pv.read(filename)
    idx = closest_point_index(grid.points, target)
    u = get_array(grid, ["f_real", "displacement_real", "f"])
    if u.ndim == 1:
        return float(abs(u[idx]))
    return float(np.linalg.norm(u[idx]))


def pressure_abs_at_file(filename, target):
    grid = pv.read(filename)
    idx = closest_point_index(grid.points, target)
    p = get_array(grid, ["acoustic_pressure_real", "f_real", "acoustic_pressure", "f"])
    return float(abs(p[idx]))


def main():
    parser = argparse.ArgumentParser(
        description="Extract pressure and displacement responses at anatomical target points."
    )
    parser.add_argument(
        "--coordinates",
        default="target_coordinates.csv",
        help="CSV file containing target coordinates.",
    )
    parser.add_argument(
        "--targets",
        nargs="+",
        default=["columella", "footplate", "prootic"],
        help="Target point names to extract from the coordinate CSV.",
    )
    parser.add_argument("--f-min", type=float, default=10.0, help="Minimum frequency in Hz.")
    parser.add_argument("--f-max", type=float, default=1000.0, help="Maximum frequency in Hz.")
    parser.add_argument("--f-step", type=float, default=50.0, help="Frequency step in Hz.")
    parser.add_argument(
        "--output-csv",
        default="target_response.csv",
        help="Output CSV file.",
    )
    parser.add_argument(
        "--output-plot",
        default="target_response_plot.png",
        help="Output plot file.",
    )
    args = parser.parse_args()

    targets = read_points_from_csv(args.coordinates, args.targets)
    frequencies = np.arange(args.f_min, args.f_max + 1, args.f_step)

    rows = []
    for idx, freq in enumerate(frequencies):
        pressure_file = Path(f"pressure_p0_{idx:06d}.vtu")
        displacement_file = Path(f"displacement_p0_{idx:06d}.vtu")

        if not pressure_file.exists() or not displacement_file.exists():
            print(f"Skipping {freq:g} Hz: missing {pressure_file} or {displacement_file}")
            continue

        for point_name, target in targets.items():
            pressure = pressure_abs_at_file(pressure_file, target)
            displacement = displacement_norm_at_file(displacement_file, target)
            rows.append(
                {
                    "frequency_hz": freq,
                    "point_name": point_name,
                    "pressure_abs": pressure,
                    "displacement_abs": displacement,
                }
            )
            print(
                f"{freq:g} Hz, {point_name}: "
                f"pressure={pressure:.6e}, displacement={displacement:.6e}"
            )

    with open(args.output_csv, "w", newline="") as f:
        fieldnames = ["frequency_hz", "point_name", "pressure_abs", "displacement_abs"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    if rows:
        fig, ax = plt.subplots(figsize=(7, 5))
        for point_name in args.targets:
            xs = [r["frequency_hz"] for r in rows if r["point_name"] == point_name]
            ys = [r["displacement_abs"] for r in rows if r["point_name"] == point_name]
            ax.plot(xs, ys, marker="o", label=point_name)

        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel("Displacement amplitude")
        ax.set_title("Target-point displacement response")
        ax.legend()
        fig.tight_layout()
        fig.savefig(args.output_plot, dpi=300)
        plt.close(fig)

    print(f"Wrote {args.output_csv}")
    if rows:
        print(f"Wrote {args.output_plot}")


if __name__ == "__main__":
    main()
