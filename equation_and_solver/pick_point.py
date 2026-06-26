#!/usr/bin/env python3
"""
Pick a point coordinate from a mesh.

This utility opens an interactive PyVista window and prints the coordinates of
the selected point. The printed coordinates can be used as constraint or output
locations in the solver and post-processing scripts.

Default input
-------------
output.msh
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pyvista as pv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Interactively pick a point from a mesh.")
    parser.add_argument("--mesh", default="output.msh", help="Input mesh file.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not Path(args.mesh).exists():
        raise FileNotFoundError(f"Mesh file not found: {args.mesh}")

    mesh = pv.read(args.mesh)

    def on_pick(point, picker) -> None:  # PyVista callback signature
        print("Picked point coordinates:", point)

    plotter = pv.Plotter()
    plotter.add_mesh(mesh, show_edges=True, color="lightgray")
    plotter.enable_point_picking(
        callback=on_pick,
        use_picker=True,
        show_message=True,
    )
    plotter.show()


if __name__ == "__main__":
    main()
