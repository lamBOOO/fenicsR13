import matplotlib.pyplot as plt

# Get data from Paraview first
h = [
    0.9886573325052778,
    0.6340332990709842,
    0.32904683851469807,
    0.16754966839339377,
    0.08734460120995041
]
l2errorsquaredx = [
    0.0126144,
    0.00394172,
    0.000697761,
    0.000103676,
    2.16907e-05
]
l2errorsquaredy = [
    0.00711858,
    0.00222254,
    0.000632597,
    9.26432e-05,
    1.69494e-05
]

plt.loglog(h, l2errorsquaredx, 'o-', label='x')
plt.loglog(h, l2errorsquaredy, 'o-', label='y')
plt.xlabel("max(h)")
plt.ylabel("error measure")
plt.title("int_Omega ( div(sigma) - grad(p) )^2 dx")
plt.legend(title='Component:')
plt.show()
plt.savefig("out.pdf", dpi=150)
