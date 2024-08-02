# TODO

- Make true 2d with 1/2 instead of 1/3 in STF for example
- Allow for hard spheres molecules
- Make all print statements also written to a logfile
- Make Lagrange multiplier method like https://fenicsproject.org/olddocs/dolfin/1.6.0/python/demo/documented/neumann-poisson/python/documentation.html
- error calculation takes long, compare https://git.rwth-aachen.de/lamBOO/fenicsR13/-/pipelines/1408567 vs https://git.rwth-aachen.de/lamBOO/fenicsR13/-/pipelines/1406382
- r13 annulus tests do not converge nicely: esol correct? compatibility condition fulfilled? but bcs are strange compared to the 2d analoge. get Mathematica notebooks and redo.
- Check r13_3d cases with exact solution from Adithya, why is there no velocity at the outer interface although it was set in Mathemtica?
  - Move also to the sphere tests and do i=6 or at least i=5.5 to further confirm
