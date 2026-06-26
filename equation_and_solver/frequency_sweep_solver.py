"""
Frequency-domain finite-element solver for the manuscript:

    Origin of the snake ear

This script runs frequency-sweep simulations using FEniCSx/dolfinx.

Required input files
--------------------
output.msh
    Volumetric mesh generated with Gmsh.

selected_vibration_face_ids.txt
    Face IDs defining the vibration boundary.

absorb_face_ids.txt
    Face IDs defining the absorbing boundary.

target_coordinates.csv
    CSV file containing the fixed point used for the Dirichlet boundary
    condition. The required row is:

        fixed_condition,x,y,z

Outputs
-------
displacement_p0_*.vtu
    Complex displacement field for each simulated frequency.

pressure_p0_*.vtu
    Complex acoustic pressure field for each simulated frequency.
"""

import argparse
import csv
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import mesh, fem, io
from dolfinx.io import gmshio
from dolfinx.fem.petsc import LinearProblem
import ufl


def read_point_from_csv(csv_path, point_name):
    """Read a named 3-D point from target_coordinates.csv."""
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Coordinate file not found: {csv_path}")

    with csv_path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("point_name") == point_name:
                return float(row["x"]), float(row["y"]), float(row["z"])

    raise ValueError(f"Point '{point_name}' was not found in {csv_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Run frequency-domain finite-element acoustic-response simulations."
    )
    parser.add_argument("--mesh", default="output.msh", help="Input Gmsh mesh file.")
    parser.add_argument(
        "--coordinates",
        default="target_coordinates.csv",
        help="CSV file containing fixed and measurement coordinates.",
    )
    parser.add_argument(
        "--fixed-point-name",
        default="fixed_condition",
        help="Point name in the coordinate CSV used for the Dirichlet condition.",
    )
    parser.add_argument(
        "--vibration-faces",
        default="selected_vibration_face_ids.txt",
        help="Text file containing vibration-boundary face IDs.",
    )
    parser.add_argument(
        "--absorb-faces",
        default="absorb_face_ids.txt",
        help="Text file containing absorbing-boundary face IDs.",
    )
    parser.add_argument("--f-min", type=float, default=10.0, help="Minimum frequency in Hz.")
    parser.add_argument("--f-max", type=float, default=1000.0, help="Maximum frequency in Hz.")
    parser.add_argument("--f-step", type=float, default=50.0, help="Frequency step in Hz.")
    parser.add_argument(
        "--pressure",
        type=float,
        default=63.2,
        help="Applied pressure amplitude on the vibration boundary.",
    )
    args = parser.parse_args()

    comm = MPI.COMM_WORLD

    domain, cell_tags, facet_tags = gmshio.read_from_msh(args.mesh, comm, rank=0, gdim=3)
    tdim = domain.topology.dim
    fdim = tdim - 1
    domain.topology.create_connectivity(fdim, tdim)
    domain.topology.create_connectivity(fdim, 0)

    # Material parameters
    E = 31.3e9
    nu = 0.3
    rho = 2400.0
    mu = E / (2 * (1 + nu))
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
    K_bulk = E / (3 * (1 - 2 * nu))
    c = 3512.0
    Z0 = rho * c
    alpha = 0.93
    R_abs = np.sqrt(1 - alpha)
    Z_r = (1 + R_abs) / (1 - R_abs)
    Z_a = Z_r * Z0

    V = fem.VectorFunctionSpace(domain, ("Lagrange", 1))
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    # Fixed point for the Dirichlet condition
    X0, Y0, Z0 = read_point_from_csv(args.coordinates, args.fixed_point_name)

    def point_selector(x):
        return (
            np.isclose(x[0], X0, atol=1e-3)
            & np.isclose(x[1], Y0, atol=1e-3)
            & np.isclose(x[2], Z0, atol=1e-3)
        )

    fixed_points = mesh.locate_entities(domain, 0, point_selector)
    if len(fixed_points) == 0:
        raise RuntimeError(
            f"No mesh vertex found near fixed point '{args.fixed_point_name}' "
            f"at ({X0}, {Y0}, {Z0})."
        )

    dofs = fem.locate_dofs_topological(V, 0, fixed_points)
    bc_value = fem.Function(V)
    bc_value.x.array[:] = 0.0 + 0.0j
    bc = fem.dirichletbc(bc_value, dofs)

    # Vibration boundary
    vibration_face_ids = np.loadtxt(args.vibration_faces, dtype=np.int32, ndmin=1)
    vibration_facets = []
    for f in range(domain.topology.index_map(fdim).size_local):
        linked_cells = domain.topology.connectivity(fdim, tdim).links(f)
        if any(cell in vibration_face_ids for cell in linked_cells):
            vibration_facets.append(f)

    vibration_facets = np.array(vibration_facets, dtype=np.int32)
    vibration_tag_value = 77
    vibration_tags = np.full(vibration_facets.shape, vibration_tag_value, dtype=np.int32)
    vibration_meshtags = mesh.meshtags(domain, fdim, vibration_facets, vibration_tags)
    ds = ufl.Measure("ds", domain=domain, subdomain_data=vibration_meshtags)

    # Absorbing boundary
    absorb_ids = np.loadtxt(args.absorb_faces, dtype=np.int32, ndmin=1)
    absorb_facets = []
    for f in range(domain.topology.index_map(fdim).size_local):
        linked_cells = domain.topology.connectivity(fdim, tdim).links(f)
        if any(cell in absorb_ids for cell in linked_cells):
            absorb_facets.append(f)

    absorb_facets = np.array(absorb_facets, dtype=np.int32)
    absorb_tags = np.full(absorb_facets.shape, 99, dtype=np.int32)
    absorb_meshtags = mesh.meshtags(domain, fdim, absorb_facets, absorb_tags)
    ds_absorb = ufl.Measure("ds", domain=domain, subdomain_data=absorb_meshtags)

    frequencies = np.arange(args.f_min, args.f_max + 1, args.f_step)

    for idx, freq in enumerate(frequencies):
        omega = 2 * np.pi * freq
        gamma = fem.Constant(domain, PETSc.ScalarType(1j * omega * rho / Z_a))

        n = ufl.FacetNormal(domain)
        p0 = fem.Constant(domain, PETSc.ScalarType(args.pressure))
        g = -p0 * n
        L = ufl.inner(g, v) * ds(vibration_tag_value)

        a_absorb = gamma * ufl.inner(u, v) * ds_absorb(99)
        a = (
            ufl.inner(
                lmbda * ufl.div(u) * ufl.Identity(3)
                + 2 * mu * ufl.sym(ufl.grad(u)),
                ufl.sym(ufl.grad(v)),
            )
            * ufl.dx
            - omega**2 * rho * ufl.inner(u, v) * ufl.dx
            + a_absorb
        )

        problem = LinearProblem(
            a,
            L,
            bcs=[bc],
            petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
        )
        uh = problem.solve()

        V_p = fem.FunctionSpace(domain, ("Lagrange", 1))
        p = fem.Function(V_p, name="acoustic_pressure")
        q = ufl.TestFunction(V_p)
        divu = ufl.div(uh)

        a_proj = ufl.inner(ufl.TrialFunction(V_p), q) * ufl.dx
        L_proj = ufl.inner(-K_bulk * divu, q) * ufl.dx

        problem_proj = fem.petsc.LinearProblem(
            a_proj,
            L_proj,
            u=p,
            petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
        )
        p_h = problem_proj.solve()

        fname_disp = f"displacement_p0_{idx:06d}.vtu"
        fname_pres = f"pressure_p0_{idx:06d}.vtu"

        with io.VTKFile(comm, fname_disp, "w") as vtkfile:
            vtkfile.write_function(uh)
        with io.VTKFile(comm, fname_pres, "w") as vtkfile:
            vtkfile.write_function(p_h)

        if comm.rank == 0:
            print(f"Completed {freq:g} Hz -> {fname_disp}, {fname_pres}")


if __name__ == "__main__":
    main()
