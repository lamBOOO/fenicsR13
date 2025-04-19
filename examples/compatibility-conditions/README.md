# Compatibility conditions

## The case of $\epsilon^w = 0$ (only relevant then)

- Then, $h(p,q) =0$ always and RHS $l_5(q) =0$ must hold.
$$
l_5(q) = \int_\Omega \dot{m} q \,\textnormal{d} \boldsymbol{x} - \int_\Gamma \left( u_n^{\textnormal{w}} - \epsilon^{\textnormal{w}} \tilde{\chi}  p^{\textnormal{w}} \right) q \,\textnormal{d} l
\overset{\mathclap{\epsilon^w = 0}}{=}
\int_\Omega \dot{m} q \,\textnormal{d} \boldsymbol{x} - \int_\Gamma u_n^{\textnormal{w}} q \,\textnormal{d} l
$$
Two cases seem problematic, $q=1$ kind of is important
1. $\int_\Omega \dot{m} 1 \,\textnormal{d} \boldsymbol{x} = 0$ works, $\int_\Omega \dot{m} 1 \,\textnormal{d} \boldsymbol{x} \neq 0$ fails (i.e. solution of $\epsilon^w=0$ and $\epsilon^w=10^{-8}$ is not close)
1. $\int_\Gamma u_n^{\textnormal{w}} 1 \,\textnormal{d} = 0$ works, $\int_\Gamma u_n^{\textnormal{w}} 1 \,\textnormal{d} \neq 0$ fails (i.e. solution of $\epsilon^w=0$ and $\epsilon^w=10^{-8}$ is not close)


## Questions
- Are the more compatibility conditions? For other equations, or other variables?
- Why is only $q=1$ relevant? $l_5$ should be zero for all test functions actually...
- Is there a connection between mass source and normal velocities? (241213 I remember yes)
  - I.e. even for $\dot{m} \neq0$, the outflow could be chosen compatible?
- Check all old notes...

