# Generation of test matrices to use with `HPDDM`

Run
```
bash run.sh
```

## Output

There will be several folders "out_*" containing:
- `A.binary`: result of `-ksp_view_mat binary:out_0/A.binary`
- `b.binary`: result of `-ksp_view_rhs binary:out_0/b.binary`
- `dofmap_0_i.mat`: which local DOF corresponds to a field on rank i
  - Note about the fields:
    - 0: scalar, P1 Lagrange
    - 1: 2d vector, P2 Lagrange
    - 2: scalar, P1 Lagrange
    - 3: 2d vector, P1 Lagrange
    - 4: 2d symmetric tensor (3 DOFs), P2 Lagrange
- `ownership_0_i.mat`: result of `ksp().getOperators()[0].getOwnershipRange()`
