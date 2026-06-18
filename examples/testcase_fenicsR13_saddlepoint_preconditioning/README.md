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
Riesz-map preconditioner. The pressure space is the corrected one from the
theory paper,

```text
V = H1(sigma) x H1(s) x H1~(p), Q = L2(u) x L2(theta),
```

where `H1~(p)` is mean-free in the pressure-nullspace case and otherwise
equivalent to `H1(p)` through the wall-pressure boundary control. This corrects
the draft note's `L2(p)` typo.

The strongest implementation keeps the pressure and stress in one `p_sigma`
block so that the equivalent Riesz norm retains the full total-pressure
boundary term

```text
epsilon_w * chi_tilde * (p + sigma_nn) * (q + psi_nn)
```

A naive separate `p`/`sigma` fieldsplit silently removes this part from the
preconditioner. A split can still be meaningful, but then the split norm has to
be a deliberate spectrally equivalent replacement, not an accidental dropped
term.

With exact LU/MUMPS block solves, the corrected Riesz run gives:

| mesh  | MINRES iterations |
|-------|------------------:|
| ring1 |               175 |
| ring2 |               228 |
| ring3 |               268 |
| ring4 |               236 |
| ring5 |               209 |

So the direct-subblock Riesz preconditioner does plateau. This is the behavior
predicted by Mardal-Winther operator preconditioning for a stable
discretization: use the Riesz map of the Hilbert norms dictated by
well-posedness, not a Schur approximation. The remaining production problem is
to replace the exact subblock inverses by spectrally equivalent scalable
solvers.

## Split `p`/`sigma` Variant

`input/r13_riesz_minres_lu_split.yml` tests the Mardal-Winther idea one step
further: replace the coupled `p_sigma` Riesz block by a product Riesz map. It
keeps boundary control by diagonalizing the total-pressure term as

```text
epsilon_w * chi_tilde * (p*q + sigma_nn*psi_nn)
```

and then solves `p` and `sigma` as independent fields. This is not the exact
same Riesz map, and the equivalence constants are worse, especially for large
wall weights where the coupled term naturally allows the cancellation
`p ~= -sigma_nn`. But the constants are still independent of mesh size in this
test:

| mesh  | MINRES iterations |
|-------|------------------:|
| ring1 |               536 |
| ring2 |               673 |
| ring3 |               632 |
| ring4 |               549 |

So yes: splitting the `p_sigma` block appears to retain mesh robustness here.
The tradeoff is a much larger iteration constant than the coupled `p_sigma`
Riesz block. That makes the split version interesting for scalable AMG
experiments, but the coupled version remains the cleaner ideal block solve.

## Iterative Block Solves

`input/r13_riesz_minres_iterative_split.yml` keeps MINRES and replaces the exact
split Riesz subblock solves by PETSc preconditioners:

- Jacobi for the L2 mass blocks `theta` and `u`
- GAMG for the H1-type blocks `s`, `p`, and `sigma`

This follows the older working `r13_riesz_minres.yml` style, but uses the
corrected `theory_split` Riesz norm instead of the legacy diagonal pressure
choice. Docker run on ring1-ring4:

| mesh  | MINRES iterations |
|-------|------------------:|
| ring1 |               679 |
| ring2 |              1005 |
| ring3 |              1095 |
| ring4 |              1146 |

The result is still qualitatively mesh-flat, but the constants are worse than
the exact split LU block solves. So the operator-preconditioning story is intact:
iterative block solvers preserve mesh robustness only to the extent that their
subblock preconditioners are spectrally equivalent with acceptable constants.
Here GAMG is usable, but not yet a sharp block inverse.

## Schur Diagnostic

`input/r13_schur_full_lu.yml` is only an exact block-factorization diagnostic. It
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

This confirms that the assembled matrix and parity flip are algebraically
consistent, but it is intentionally expensive and is not the operator
preconditioner advocated by the theory.

Do not confuse this diagnostic split with directly pivoting on the continuous
paper block `V=(sigma,s,p)`: `A` is only coercive on `ker B`, so that block is
not a robust standalone Schur pivot in the discrete matrix.

## Practical Direction

For a scalable preconditioner, the target is to approximate the Riesz blocks
well enough: multigrid for the H1-type `s` block, either coupled `p_sigma`
multigrid or the split product-norm variant for pressure/stress, and
lumped/Jacobi-quality solves for the L2 multiplier blocks. For MINRES, the
approximate block preconditioner must remain SPD; the split GAMG/Jacobi version
is the current practical MINRES candidate. The Schur file stays as a sanity
check, not as the proposed mesh-independent method.
