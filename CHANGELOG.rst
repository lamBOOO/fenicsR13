Change log
----------

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
