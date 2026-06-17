# R13 Saddle-Point Preconditioning

## Problem

The original Schur/`selfp` setup was not mesh independent:

| mesh  | outer iterations |
|-------|-----------------:|
| ring1 |              421 |
| ring2 |              868 |

Direct solves inside the visible subblocks are not enough if the
preconditioner drops the coupling that controls the discrete saddle system.

## Theory Riesz Preconditioner

`input/r13_riesz_minres_lu.yml` uses the parity-symmetrized system and an SPD
Riesz-map preconditioner. The implementation follows the paper grouping

```text
V = (sigma, s, p), Q = (u, theta)
```

and fixes the important pressure/stress issue: the pressure and stress must be
solved as one `p_sigma` block so that the boundary term

```text
epsilon_w * chi_tilde * (p + sigma_nn) * (q + psi_nn)
```

is retained. Splitting `p` and `sigma` separately silently removes this
total-pressure part from the preconditioner.

With exact LU/MUMPS block solves, the corrected Riesz run gives:

| mesh  | MINRES iterations |
|-------|------------------:|
| ring1 |               175 |
| ring2 |               228 |
| ring3 |               268 |

This is much better than the old exact-block Riesz counts, but it is not a
literal one-iteration algebraic ideal. The remaining drift is consistent with
the chosen discrete pair not giving the exact Schur action
`B R_V^{-1} B^T`.

## Algebraic Ideal

`input/r13_schur_full_lu.yml` is the exact block-factorization diagnostic. It
uses:

- `solver_options.symmetrize: True`
- PETSc `pc_fieldsplit_detect_saddle_point`
- PETSc `pc_fieldsplit_type: schur`
- PETSc `pc_fieldsplit_schur_precondition: full`
- LU/MUMPS on both PETSc-detected subblocks

Docker run on ring1-ring3:

| mesh  | GMRES iterations |
|-------|-----------------:|
| ring1 |                1 |
| ring2 |                1 |
| ring3 |                1 |

This confirms that the truly mesh-independent algebraic ideal must include the
full discrete Schur complement. It is intentionally expensive and is a
diagnostic/reference preconditioner, not a scalable production method.

Do not confuse this diagnostic split with directly pivoting on the continuous
paper block `V=(sigma,s,p)`: `A` is only coercive on `ker B`, so that block is
not a robust standalone Schur pivot in the discrete matrix.

## Practical Direction

For a scalable preconditioner, the target is now clear: approximate the
discrete Schur action `B R_V^{-1} B^T` well enough, or use a finite element
pair/stabilization with a mesh-uniform discrete inf-sup constant. Improving
the direct solves inside independent Riesz subblocks cannot fix the missing
Schur coupling by itself.
