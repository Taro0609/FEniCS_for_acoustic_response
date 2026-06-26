#!/usr/bin/env python3
"""
Select absorbing boundary faces from a Gmsh mesh.

This script opens an interactive PyVista window, allows repeated cell picking
on the extracted surface, and writes the selected face IDs to a text file.

Default input
-------------
output.msh

Default outputs
---------------
absorb_face_ids.txt
absorb_face.vtu
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Interactively select absorbing boundary faces from a mesh."
    )
    parser.add_argument("--mesh", default="output.msh", help="Input mesh file.")
    parser.add_argument(
        "--ids-out",
        default="absorb_face_ids.txt",
        help="Output text file containing selected face IDs.",
    )
    parser.add_argument(
        "--surface-out",
        default="absorb_face.vtu",
        help="Output VTU file containing the selected surface cells.",
    )
    return parser.parse_args()


def pick_faces(mesh_file: str) -> set[int]:
    mesh = pv.read(mesh_file)
    surface = mesh.extract_surface()
    picked_faces: set[int] = set()

    print("Interactive face selection")
    print("Use cell picking in the PyVista window. Close the window to continue.")
    print("Repeat the selection when prompted, or enter 'n' to finish.")

    while True:
        plotter = pv.Plotter()
        picked_this_round: list[bool] = []

        def callback(picked: pv.DataSet) -> None:
            if picked.n_cells == 0:
                print("No cells selected.")
                return

            if "vtkOriginalCellIds" in picked.array_names:
                ids = picked["vtkOriginalCellIds"]
            else:
                ids = np.arange(picked.n_cells)

            picked_faces.update(int(i) for i in ids)
            picked_this_round.append(True)
            print(f"Selected {len(ids)} cells in this round; total = {len(picked_faces)}.")

        plotter.add_mesh(surface, show_edges=True, color="lightgray")
        plotter.enable_cell_picking(
            callback=callback,
            through=True,
            show=True,
            style="surface",
        )
        plotter.show()

        if not picked_this_round:
            print("No cells were selected in this round.")

        answer = input("Select more faces? [y/N]: ").strip().lower()
        if answer != "y":
            break

    return picked_faces


def main() -> None:
    args = parse_args()

    if not Path(args.mesh).exists():
        raise FileNotFoundError(f"Mesh file not found: {args.mesh}")

    picked_faces = pick_faces(args.mesh)

    if not picked_faces:
        print("No faces selected. No output files written.")
        return

    mesh = pv.read(args.mesh)
    surface = mesh.extract_surface()
    selected = surface.extract_cells(sorted(picked_faces))
    selected.save(args.surface_out)

    np.savetxt(args.ids_out, np.array(sorted(picked_faces), dtype=np.int32), fmt="%d")

    print(f"Wrote selected face IDs to: {args.ids_out}")
    print(f"Wrote selected surface cells to: {args.surface_out}")


if __name__ == "__main__":
    main()
