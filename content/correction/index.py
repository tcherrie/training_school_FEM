# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.6
#   kernelspec:
#     display_name: ngsolve_latest
#     language: python
#     name: python3
# ---

# %% [markdown]
# # FEM with NGSolve - Tutorial Index
#
# ## NUMELEC 2026 Training School (2 days)
# ---
#
# ## Day 1 - Foundations: 2D magnetostatics
#
# | Notebook | Description |
# |---|---|
# | [**Day 1 A - First steps in FEM universe**](day1_a_first_steps.ipynb) | Mesh, function spaces, variational formulations, assembly, solving |
# | [**Day 1 B - 2D Magnetostatics**](day1_b_2Dmagnetostatics.ipynb) | Magnetostatic problems, nonlinear materials, Newton solver |
# | [**Day 1 C - Mesh refinement**](day1_c_mesh_refinement.ipynb) | Adaptive mesh refinement, error estimators, convergence studies |
#
# ---
#
# ## Day 2 - Advanced Topics: 3D magnetoquasistatics
#
# | Notebook | Description |
# |---|---|
# | [**Day 2 A - Direct and iterative solvers**](day2_a_direct_and_iterative_solvers.ipynb) | Direct (Cholesky, SuperLU) vs iterative (CG, GMRES, QMR, MINRES) solvers |
# | [**Day 2 B - Gauge methods**](day2_b_gauge.ipynb) | Singularity handling in H(curl) spaces: mass term and Coulomb gauge |
# | [**Day 2 C - 3D eddy currents**](day2_c_3D_eddy_current.ipynb) | Eddy current formulation, TEAM benchmark problem 15 |
#
# ---
#
# ## Utility modules
#
# The following modules live in the `utils/` folder and are imported by the notebooks:
#
# | Module | Contents |
# |---|---|
# | `utils.linalg` | Sparse matrix helpers: reordering, spectrum, sparsity visualisation |
# | `utils.myGeometries` | 2D and 3D geometry builders (square, disk, inductor, TEAM 15...) |
# | `utils.mySolvers` | `mySolver` (linear) and `solveNL` (nonlinear Newton with line search) |
# | `utils.ngcotree` | Co-tree basis construction for H(curl) problems |
#
# ---
#
# ## Quick navigation
#
# - [report.md](report.md) - Full code review, per-file status, and issue list
# - [Day 1](#day-1--foundations) to Start here for the basics
# - [Day 2](#day-2--advanced-topics) to Dive into advanced numerical methods
#
# > **Tip:** Open any notebook directly by clicking its link above, or run it from the command line:
# > ```
# > python day1_a_first_steps.py
# > python day2_c_3D_eddy_current.py
# > ```

# %% [markdown]
# ### Authors
# - Théodore Cherrière (CentraleSupélec, GeePs)
# - Innocent Niyonzima (Université Grenoble Alpes, G2eLab)
# - Zuqi Tang (Université de Lille, L2EP)
# - Guillaume Wasselynck (Nantes Université, IREENA)

# %% [markdown]
# ### Note
#
# This tutorial relies on [NGSolve](https://ngsolve.org/), along with [JupyterLite](https://jupyterlite.readthedocs.io/en/stable/). The jupyterlite implementation is based on the this [template GitHub repository](https://github.com/NGSolve/jupyterlite_ngsolve).
#
# ### References
#
# - J. Schöberl (TU Wien), **The Finite Element Method (iFEM), interactive NGSolve notebooks**: [https://jschoeberl.github.io/iFEM/](https://jschoeberl.github.io/iFEM/)
#
# - [TEAM Workshop Problem 15 (1997)](https://www.compumag.org/jsite/images/stories/TEAM/problem15.pdf)
#
