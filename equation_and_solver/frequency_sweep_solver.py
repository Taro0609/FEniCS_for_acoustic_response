"""
Frequency-domain finite-element solver for the manuscript:

    Origin of the snake ear

This script performs frequency-sweep simulations of acoustic responses
using FEniCSx/dolfinx.

Input files
-----------
output.msh
    Volumetric mesh generated from Gmsh.

selected_vibration_face_ids.txt
    Surface IDs defining the vibration boundary.

absorb_face_ids.txt
    Surface IDs defining the absorbing boundary.

Output
------
displacement_p0_*.vtu
    Complex displacement field for each simulated frequency.

pressure_p0_*.vtu
    Complex acoustic pressure field for each simulated frequency.
"""

import os
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import mesh, fem, io
from dolfinx.io import gmshio
from dolfinx.fem.petsc import LinearProblem
import ufl

#1_mesh_loading
comm = MPI.COMM_WORLD
domain, cell_tags, facet_tags = gmshio.read_from_msh("output.msh", comm, rank=0, gdim=3)
tdim = domain.topology.dim
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)
domain.topology.create_connectivity(fdim, 0)

#2_material_properties
E = 31.3e9
nu = 0.3
rho = 2400.0
mu = E / (2 * (1 + nu))
lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
K_bulk = E / (3 * (1 - 2 * nu))
c = 3512  # sound velocity
Z0 = rho * c
alpha = 0.93
R_abs = np.sqrt(1 - alpha)
Z_r = (1 + R_abs) / (1 - R_abs)
Z_a = Z_r * Z0

#3_function_space
V = fem.VectorFunctionSpace(domain, ("Lagrange", 1))
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

#4_fixed_condition
X0, Y0, Z0 = 4.698681, 18.07331448, 45.74733412  #ventral tip of the prootic

def point_selector(x):
    return (
        np.isclose(x[0], X0, atol=1e-3) &
        np.isclose(x[1], Y0, atol=1e-3) &
        np.isclose(x[2], Z0, atol=1e-3)
    )

fixed_points = mesh.locate_entities(domain, 0, point_selector)
dofs = fem.locate_dofs_topological(V, 0, fixed_points)
bc_value = fem.Function(V)
bc_value.x.array[:] = 0.0 + 0.0j
bc = fem.dirichletbc(bc_value, dofs)

#5_Neumann_condition
vibration_face_ids = np.loadtxt("selected_vibration_face_ids.txt", dtype=np.int32)
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

#6_absorption_condition
absorb_ids = np.loadtxt("absorb_face_ids.txt", dtype=np.int32, ndmin=1)
absorb_facets = []
for f in range(domain.topology.index_map(fdim).size_local):
    linked_cells = domain.topology.connectivity(fdim, tdim).links(f)
    if any(cell in absorb_ids for cell in linked_cells):
        absorb_facets.append(f)
absorb_facets = np.array(absorb_facets, dtype=np.int32)
absorb_tags = np.full(absorb_facets.shape, 99, dtype=np.int32)
absorb_meshtags = mesh.meshtags(domain, fdim, absorb_facets, absorb_tags)
ds_absorb = ufl.Measure("ds", domain=domain, subdomain_data=absorb_meshtags)

#7_frequency_sweep_loop
frequencies = np.arange(10, 1000 + 1, 50)
for idx, f in enumerate(frequencies):
 omega = 2 * np.pi * f
 gamma = fem.Constant(domain, PETSc.ScalarType(1j * omega * rho / Z_a))
 n = ufl.FacetNormal(domain)
 p0 = fem.Constant(domain, PETSc.ScalarType(63.2))  #130dB
 g = -p0 * n
 L = ufl.inner(g, v) * ds(vibration_tag_value)
 a_absorb = gamma * ufl.inner(u, v) * ds_absorb(99)
 a = (
    		ufl.inner(lmbda*ufl.div(u)*ufl.Identity(3) + 2*mu*ufl.sym(ufl.grad(u)), ufl.sym(ufl.grad(v))) * ufl.dx
    - omega**2 * rho * ufl.inner(u, v) * ufl.dx + a_absorb
		)

 problem = LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})

 uh = problem.solve()
 V_p = fem.FunctionSpace(domain, ("Lagrange", 1))
 p = fem.Function(V_p, name="acoustic_pressure")
 q = ufl.TestFunction(V_p)
 divu = ufl.div(uh)
 a_proj = ufl.inner(ufl.TrialFunction(V_p), q) * ufl.dx
 L_proj = ufl.inner(-K_bulk * divu, q) * ufl.dx
 problem_proj = fem.petsc.LinearProblem(a_proj, L_proj, u=p,petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
 p_h = problem_proj.solve()
 fname_disp = f"displacement_p0_{idx:06d}.vtu"
 fname_pres = f"pressure_p0_{idx:06d}.vtu"

 with io.VTKFile(comm, fname_disp, "w") as vtkfile:
  vtkfile.write_function(uh)
 with io.VTKFile(comm, fname_pres, "w") as vtkfile:
  vtkfile.write_function(p_h)

 if comm.rank == 0:
  print(f"Completed {f} Hz → {fname_disp}, {fname_pres}")

