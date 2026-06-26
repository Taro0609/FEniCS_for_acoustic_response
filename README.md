Origin of the Snake Ear: Finite Element Analysis

This repository contains the scripts used for the finite element analyses described in the manuscript:

Origin of the snake ear

The code reproduces the acoustic response simulations used in the study. Statistical analyses and morphometric measurements are not included in this repository unless otherwise stated.

Repository contents
mesh/
input/
output/
scripts/
mesh/ : finite-element meshes
input/ : input files
output/ : simulation results
scripts/ : Python scripts for running the analyses

Requirements

The simulations were developed using:

Python 3
FEniCSx (dolfinx)
PETSc
Gmsh
PyVista
Docker

The analyses were run under Docker to ensure a consistent computational environment.

Installation

Install Docker and clone this repository.

The Docker image provides all required software dependencies.

Running the simulations

Run the Python scripts in the scripts/ directory after placing the required mesh and input files in the appropriate folders.

Simulation outputs are written to the output/ directory.

Reproducing the manuscript

The repository contains the scripts used to generate the finite element simulations reported in the manuscript.

The meshes and input files correspond to the specimens analysed in the study.

License

MIT License
