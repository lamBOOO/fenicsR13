# Readme

## Problem
- no mesh independence (doubling of its. if h |-> h/2):
mesh 1: 421 outer its.
mesh 2: 868 outer its.

## PETSc options (see input file)
```
petsc_options:
  ksp_type: tfqmr
  ksp_rtol: 1E-7
  ksp_max_it: 100000
  ksp_view:
  ksp_monitor_true_residual:
  pc_type: fieldsplit
  pc_fieldsplit_detect_saddle_point:
  # option 1:
  pc_fieldsplit_type: schur
  pc_fieldsplit_schur_fact_type: full
  pc_fieldsplit_schur_precondition: selfp
  fieldsplit_0_ksp_type: preonly
  fieldsplit_0_pc_type: gamg
  fieldsplit_0_pc_gamg_type: classical
  fieldsplit_1_ksp_type: preonly
  fieldsplit_1_pc_type: jacobi
```

## Matrices and Vectors
Can be found in `matrix_and_vector` folder. A_i,b_i correspond to the i-th mesh, A_0,b_0 to the coarsest mesh.

## Paraview output
is in `output` folder. Might not be relevant.

## Internal notes
- Test case adapted from https://git.rwth-aachen.de/lamBOO/fenicsR13/-/blob/master/tests/2d_r13/inputs/r13_1_coeffs_nosources_norot_inflow_p1p2p1p1p2_nostab.yml?ref_type=heads
- Used Schur + MG preconditioner, enabled "write_systemmatrix" to write the system matrix to a file
