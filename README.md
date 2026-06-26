# FEniCS acoustic-response simulations

Code used for the finite-element acoustic-response simulations in the manuscript:

**Origin of the snake ear**

This repository contains the scripts used to prepare meshes, define boundary conditions, run frequency-domain simulations in FEniCSx/dolfinx, and export displacement and pressure fields.

The simulations were used to estimate frequency-dependent mechanical responses of the jaw-ear region under jaw-borne loading.

## Repository contents

* `meshing/`
  Scripts used for STL preprocessing and generation of volumetric meshes with Gmsh.

* `helmholtz_equation/`
  Scripts related to the acoustic Helmholtz component of the model.

* `equation_and_solver/`
  Scripts for material-property definition, boundary conditions, frequency-domain solving, and output export.

## System requirements

The scripts were run using:

* Python 3
* Docker
* FEniCSx / dolfinx v0.7.3
* Gmsh v4.11
* PETSc
* PyVista
* NumPy
* mpi4py
* petsc4py

The simulations were run in a Docker container using the complex-mode dolfinx environment.

## Installation

Install Docker and clone this repository.

```bash
git clone https://github.com/Taro0609/FEniCS_for_acoustic_response.git
cd FEniCS_for_acoustic_response
```

Start the dolfinx container:

```bash
docker run -ti -v ${PWD}:/home/fenics/shared -w /home/fenics/shared -e DISPLAY=host.docker.internal:0.0 --network=host dolfinx/dolfinx:v0.7.3 bash -c "source dolfinx-complex-mode && bash"
```

The installation time depends mainly on the Docker image download. Once the image is available locally, the container starts within a few minutes on a standard desktop computer.

## Workflow

The complete workflow consists of five steps.

1. Generate a volumetric mesh (`output.msh`) from the segmented STL model.
2. Select the vibration boundary using `pick_vibration_face.py`.
3. Select the absorbing boundary using `pick_absorption_face.py`.
4. Run `frequency_sweep_solver.py`.
5. Analyse simulation outputs using `target_plot.py`.

## Required input files

The solver requires the following files in the working directory:

- `output.msh`  
  Volumetric mesh generated with Gmsh.

- `selected_vibration_face_ids.txt`  
  Face IDs defining the vibration boundary.

- `absorb_face_ids.txt`  
  Face IDs defining the absorbing boundary.

## Running the simulations

Place the required mesh and boundary-condition files in the working directory.

Then run the solver script from inside the Docker container:

```bash
python3 frequency_sweep_solver.py
```

The solver runs a frequency sweep and writes displacement and pressure outputs as `.vtu` files.

## Expected outputs

The main outputs are:

* `displacement_p0_*.vtu`
* `pressure_p0_*.vtu`

These files contain complex-valued displacement and pressure fields for each simulated frequency.

The displacement outputs were used to quantify frequency-dependent displacement amplitudes at anatomical output sites, including the columella, footplate and prootic.

## Notes

Specimen-specific meshes and manually selected boundary-condition files are required to reproduce the full analysis.

Large CT-derived mesh files are not stored directly in this repository unless otherwise specified.

## License

MIT License

## Citation

If you use this code, please cite:

Nojiri T. et al.
**Origin of the snake ear.**

## Contact

Taro Nojiri
