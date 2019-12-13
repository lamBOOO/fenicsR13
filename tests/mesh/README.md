# Meshes

- Two versions: Random ring, antisymmetric ring with flat upper corner of inner circle
- Random ring mesh version produces actually indefinite (one eigenvalue slightly below zero in machine precision region) for the case "r13_1_coeffs_sources_rot_noinflow_p1p1p1p1p1_stab" using mumps solver, this means a non-unique solution? Default solver ("direct") works however..
  - Actually, all systems are indefinite if rescale pressure is needed!!!! Sometimes we have -10E-15, sometimes +10E-15