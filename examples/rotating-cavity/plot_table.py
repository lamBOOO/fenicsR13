import matplotlib.pyplot as plt
import numpy as np


gui = False
if not gui:
    plt.switch_backend("agg")
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Plot curves for table_0.csv .. table_4.csv if present
for i in range(5):
    fname = f"table_{i}.csv"
    try:
        x, y = np.loadtxt(fname, delimiter=",", unpack=True)
        xstar = 0.1 * (i + 1)
        label = (
            rf"$\mathbf{{x}}_{{\mathrm{{bot}}}}/L = "
            rf"{xstar:.1f}$"
        )
        ax.plot(x, y, "-o", label=label, markersize=4)
    except Exception as e:
        print(f"WARN: Could not load {fname}: {e}")

ax.set_xscale("log")

plt.subplots_adjust(left=0.18, right=0.95, top=0.9, bottom=0.15)
# add minor grid
ax.minorticks_on()
ax.grid(which="both", linestyle="--", linewidth=0.1)
plt.title("Rotating Cavity: Kinetic Energy vs Knudsen Number")
plt.xlabel(r"Knudsen number $\mathrm{Kn}$")
plt.ylabel(
    (
        r"kinetic energy: $\frac{1}{2}\Vert \mathbf{u} \Vert^2_{L^2(\Omega)} = "
        r"\frac{1}{2} \int_{\Omega} (\mathbf{u} \cdot \mathbf{u}) d\mathbf{x}$"
    )
)
plt.legend()

gui = False
if gui:
    plt.show()
else:
    fig.savefig("fig.pdf", dpi=150)
