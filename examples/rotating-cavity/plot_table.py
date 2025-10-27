import matplotlib.pyplot as plt
import numpy as np

x, y = np.loadtxt("table.csv", delimiter=",", unpack=True)


gui = False
if not gui:
    plt.switch_backend("agg")

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(x, y, "-o", label="max ||u||_2")
ax.set_xscale("log")

plt.xlabel("Knudsen number")
plt.ylabel("max ||u||_2")
plt.title("Rotating Cavity: Max Velocity vs Knudsen Number")
plt.legend()

gui = False
if gui:
    plt.show()
else:
    fig.savefig("fig.pdf", dpi=150)
