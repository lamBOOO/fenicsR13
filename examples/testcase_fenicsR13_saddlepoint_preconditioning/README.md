# Readme

## Problem
- no mesh independence (doubling of its. if h |-> h/2):
mesh 1: 421 outer its.
mesh 2: 868 outer its.

## Solution: Riesz-map operator preconditioning (see `input/r13_riesz_minres.yml`)
Instead of Schur/`selfp`, use the well-posedness structure directly:
- `solver_options.symmetrize: True` flips the sign of the even-parity test
  rows (kappa, psi, q), turning the assembled system into a *symmetric*
  indefinite saddle point matrix (verified at runtime: `||K - K^T||_F ~ 1e-13`)
  => MINRES applicable.
- `solver_options.riesz_preconditioner: True` assembles a separate SPD
  block-diagonal preconditioner from the Riesz maps of the well-posedness
  norms (Mardal-Winther): the coercive diagonal blocks `a(s,r)`, `d(sigma,psi)`
  (H1-type, incl. Kn-weights and boundary terms) for s and sigma, L2 mass for
  theta and u, and H1 (mass + grad) for p. CIP jumps enter P with positive
  sign when enabled.
- The solver attaches per-field index sets to PETSc `fieldsplit` (named
  splits `fieldsplit_<theta|s|p|u|sigma>_*`) and composes near-nullspaces
  onto the H1 blocks for GAMG (rigid body modes for the elasticity-like
  s-block, per-component constants for the tensor-valued sigma-block).

Result (MINRES, additive fieldsplit, GAMG/Jacobi, rtol 1e-7):

| mesh  | hmax  | iterations (serial) | iterations (mpirun -np 2) |
|-------|-------|--------------------:|--------------------------:|
| ring1 | 0.634 |                 564 |                       666 |
| ring2 | 0.329 |                 871 |                       961 |
| ring3 | 0.168 |                1036 |                      1193 |
| ring4 | 0.087 |                1264 |                      1431 |

The near-nullspace construction is MPI-parallel (modes are extracted with
`VecGetSubVector` on the field IS, matching the fieldsplit submatrix
layout; Gram-Schmidt uses global PETSc dots). The ~15% serial-to-parallel
overhead is the usual parallel AMG aggregation/smoother variation.

Iteration ratios decay (1.54 -> 1.19 -> 1.22), i.e. approaching mesh
independence, vs. clean doubling (and eventual stagnation) for the original
tfqmr + Schur/`selfp` setup. With exact block solves (cholesky on the H1
blocks) the counts are 488/727/820/907 (ring1-4, ratios 1.49/1.13/1.11,
i.e. saturating), confirming the operator preconditioning theory; the
remainder is AMG block quality plus a slight residual drift (likely the
H1 p-block heuristic, see below).

Empirical findings (relevant for the theory write-up):
- An L2 mass block for p (as the V-norm of Theorem 4.9 would suggest) yields
  h-dependent counts (~doubling per refinement, even with exact block
  solves): the implemented pressure coupling `g(p,v) = (v, grad p)` is not
  L2(p) x L2(u) bounded. An H1 p-block restores (near) mesh independence.
- CIP stabilization alone does not fix this (875/2797/4312 with CIP and
  L2 p-block, exact block solves).
- Accuracy is unaffected by the parity flip: a convergence study against the
  exact solution reproduces the expected error decay.

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
