Change log
----------

1.3 (2020-08-18)
~~~~~~~~~~~~~~~~

- Add option to use Garlerkign Least Squares (GLS) stabilization
- Allow for different chi_tilde on boundaries
- Add option for multiple mesh regions:
    - Allows for different Knudsen numbers in regions
    - Add new example: Lid-driven cavity with two mesh regions
- Optimization of CI pipeline and scripts in Gitlab:
    - Assert Flake8 compliance
- Add option to output solution vectors (for model order reduction experiments)
- Improve README with R13 equation set and list of features

1.2 (2020-02-17)
~~~~~~~~~~~~~~~~

- Change system:
    - Change RHS for energy coupling consistency
    - Add delta term
    - Add body force
    - Add div(u) coupling to energy balance (also change RHS for elimination)
    - Rename deltas and chi_tilde
- Add body froce driven channel flow
- Add P2P1P1P2P1 R13 validation case
- Revert to MUMPS solver

1.1 (2019-12-15)
~~~~~~~~~~~~~~~~

- Symmetrize system
    - Scale the equations to match
    - Introduce subfunctionals
    - Subfunctionals are equal for off-diagonal entries
    - Resolve stf-terms using orthogonality
    - Make r13 as default in formulation rather than decoupled
    - Refactor CIP as separate subfunctionals
    - Add antisymm ring, not used for now but can improve convergence
- Create package
    - Move files into separate folder and add setup.py
    - CI is changed
    - Installation through "pip install ." in toplevel
- Clean repository

1.0 (2019-09-23)
~~~~~~~~~~~~~~~~

- Improve examples
    - Remove auxiliary files, clean cases, improve output plots and change some meshes
- Add WELCOME screen to Docker container with link to website and documentation
- Introduce more variables as ``dolfin.Constant()`` for parameters study w.o. compiling
- Add time measure for linear solve
- Thesis submission version

0.5 (2019-09-13)
~~~~~~~~~~~~~~~~

- Add more examples with tests and documentation:
    - Lid-Driven Cavity
    - Channel Flow with Knudsen paradox plot
    - Knudsen pump
- Add option to perform parameter studies
- Add massflow reporting option for arbitrary BCs
- Fix P1P2P4 stress test case
- Add more printing statement to program output
- Change formulation:
    - Rename gamma to epsilon in inflow model
    - Rename tau to Knudsen to have real dimensionless equations
    - Replace sym(psi) -> psi because symmetric per definition
    - Fix stf3d2 for arbitrary
- Extend documentation:
    - Extended tutorial
    - Move legacy notes to bottom of README

0.4 (2019-08-21)
~~~~~~~~~~~~~~~~

- Finish documentation
    - Also includes some ``doctests`` to test for edge cases
- Introduce develop branch to only have major version at master branch
- Add relative error calculation
- Add channel flow example (to test for Knudsen paradox)
- Fix error calculation for higher-order Ansatz function
    - The previous error was based on DOFS (P2 elements therefore differ), the new error is based on vertex values

0.3 (2019-08-11)
~~~~~~~~~~~~~~~~

- Full linear R13 now converges
- Inflow model works
- Restructuration of BC specification
- Minor improvements in plotting routine

0.2 (2019-07-29)
~~~~~~~~~~~~~~~~

- Decoupled stress system converges
- Add separated tensor operations module
    - This was needed to implement operations on synthetic 3D tensors
- Add pytests for stress
- Add new logo
- Add more Sphinx documentation
- Restructure repository

0.1 (2019-07-17)
~~~~~~~~~~~~~~~~

- Add logo
- Add Sphinx documentation
- Add pytests
- Add Gitlab CI scripts
