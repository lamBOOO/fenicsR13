import matplotlib.pyplot as plt
import numpy as np

x, y = np.loadtxt("table.csv", delimiter=",", unpack=True)


gui = False
if not gui:
    plt.switch_backend("agg")

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(x,y, "-o", label="Channel Flow")
ax.set_xscale("log")

plt.xlabel("Knudsen number")
plt.ylabel("Dimensionless Massflow")
plt.title("Knudsen Paradox")
plt.legend()

gui = False
if gui:
    plt.show()
else:
    fig.savefig("fig.pdf", dpi=150)
