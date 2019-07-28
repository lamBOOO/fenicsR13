# Workplan/Roadmap

## Project Description

The goal of the thesis is to enable simulations of gas flows (e.g. rarefied or diluted gases) in non-equilibrium. The open-source computing platform FEniCS [1] is used to solve the corresponding partial differential equations using the Finite Element Method. The full regularized 13-moment (R13) equations are first linearized and two decoupled systems of equations are considered. Numerical solutions of these systems have to be studied compared to analytical solutions that can be found in existing work in this research field [3].
Part of this work is also a step towards the numerical analysis of more complex model equations and problems describing flows in microscopic settings, i.e.:

- generation and validation of new test problems
- inclusion of the coupling terms to solve the full linear system
- implementation of the above-listed models in the three-dimensional setting
- extension of the steady models to the transient case

References

[1]: M. S. Alnaes, J. Blechta, J. Hake, A. Johansson, B. Kehlet, A. Logg, C. Richardson, J. Ring, M. E. Rognes and G. N. Wells. The FEniCS Project Version 1.5. Archive of Numerical Software, Vol. 3 (2015).
[2]: H. Struchtrup and M. Torrilhon. Regularization of Grad's 13-Moment-Equations: Derivation and Linear Analysis. Physics of Fluids 15/9 (2003).
[3]: A. Westerkamp, M. Torrilhon. Finite element methods for the linear regularized 13-moment equations describing slow rarefied gas flows. Journal of Computational Physics, Vol. 389 (2019).

## Expected Workload Distribution

- Installation of required software and training - 1 week
- Literature research and research of previous work on this topic - 1 week
- Implementation of the decoupled heat system and comparison to the analytical solution - 2 weeks
- Implementation of the decoupled stress system and comparison to the analytical solution - 3 weeks
- Implementation of the full linear system and validation - 4 weeks
- Thesis writing and documentation of the work - 7 weeks
