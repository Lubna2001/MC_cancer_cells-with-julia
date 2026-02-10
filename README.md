# Agent based lattice tumor growth model (Julia)

This repository contains a Julia agent based lattice model for tumor growth developed from scratch for this project.  
The model simulates normal and cancer cell dynamics on a two dimensional lattice and supports experiments that vary mutation timing and random seed to study stochastic effects and spatial constraints.

---

## Repository structure
├─ src/ # Julia source code
│ ├─ MC_cancer_cell_model.jl # core model (cell states and update rules)
│ ├─ Simulate.jl # simulation driver (runs model, collects outputs)
│ ├─ Plots.jl # plotting utilities (snapshots, growth curves)
│ └─ Plots_seed6.jl # main analysis script (fixed seed, varying mutation time)
│
├─ figures/ # generated figures and subfolders from analysis/testing
├─ documents/ # LaTeX sources and report material
├─ .gitignore
└─ README.md


## Usage
For a minimal runnable example, see:

scripts/HOW_TO_USE.jl
