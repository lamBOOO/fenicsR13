# gasDynamicsFEM

[![pipeline status](https://git.rwth-aachen.de/lambert.theisen/gasdynamicsfem/badges/master/pipeline.svg)](https://git.rwth-aachen.de/lambert.theisen/gasdynamicsfem/commits/master)
[![coverage report](https://git.rwth-aachen.de/lambert.theisen/gasdynamicsfem/badges/master/coverage.svg)](https://git.rwth-aachen.de/lambert.theisen/gasdynamicsfem/commits/master)

Repository for Master thesis project regarding FEM simulations for non-equilibrium gas dynamics.

## Usage with Docker

The main folder of this repository contains a `Dockerfile` defining the used environment. Here, we used the optimized and official FEniCS Docker image and include `Gmsh` and install some requirements from the `requirements.txt`. This can take a while, especially the `Gmsh`mirror can be quite slow. To avoid very long execution commands (`docker run <..> -v <volume share> <etc..>`), a `docker-compose.yml` is used to store all these parameters. `docker-compose` acts as an wrapper for the Docker execution.

The `fenics` environment (also called *service* in the `docker-compose.yml`) first has to be build and can be executed afterwards. The steps to perform then read

```
docker-compose build fenics
docker-compose run fenics
```

The whole repository is mounted as a volume under `/home/fenics/shared` in the container and should be the default folder on startup. To execute the solver, move to the case folder (e.g. `/home/fenics/shared/cases/heatSystem`) and execute the script (e.g. `python3 heat.py`). Output files should be written in that case, e.g. to the `results` folder.

It is convenient to use a Jupyter sever or a X11 forwarding.

### Interactive Jupyter Notebooks with Microsoft's Visual Studio Code

This is the most convenient solution. A tutorial will follow.
Run a file with `%run ../../src/fenicsr13.py`

### X11 Window Forwarding on OSX

See [this guide](http://joshuamccall.com/articles/docker.html) for the programs to install. Then source the `open-macos-gui-tunnel.sh` with `. open-macos-gui-tunnel`. Afterwards, start the container and run the `change-matplotbib-backend-tkagg.sh` script to set the right `matplotlib`'s output.

### X11 Window Forwarding on Windows

This has to be studied.

## macOS Native FEniCS Installation

1. Install `miniconda` from [here]([here](https://conda.io/projects/conda/en/latest/user-guide/install/macos.html))
   1. If using `zsh`, add miniconda bins to PATH: `export PATH="$HOME/miniconda3/bin:$PATH"` to `~/.zshrc`
   2. Maybe, activation has to be done with executing `<path to miniconda>/bin/activate`
   3. Optional: Create separate coda environment: `conda creafenics-env`
2. Install FEniCS using conda: `conda install -c conda-forge fenics`
   1. Optional: Install `matplobib`: `conda install -c conda-forge matplotlib`
   2. Optional: Install `meshio`: `conda install -c mrossi meshio`
   3. Optional (for linting): `conda install pylint`
   4. Install mshr with `conda install -c conda-forge mshr`
   5. Fix macOS bug in matplotbib: `mkdir -p ~/.matplotlib; echo "backend: TkAgg" > ~/.matplotlib/matplotlibrc`
   6. XCode and command line developer tools msut be installed!
   7. Optional: Install Jupyter: `conda install -c anaconda jupyter`
   8. Optional: Install documentation system: `conda install -c anaconda sphinx`
   9. Optional: `conda install -c anaconda sympy`

## Docker Usage

TODO

### Docker Installing

1. Install Docker for OS, e.g. from: [here](https://hub.docker.com/editions/community/docker-ce-desktop-mac)

### FEniCS Install Using Docker [(see here)](https://fenics.readthedocs.io/projects/containers/en/latest/quickstart.html)

1. `curl -s https://get.fenicsproject.org | bash`
2. Run with `fenicsproject run`

### Tips/Bugs

- Matplotbib fails when having wrong backend on macOS
  - Fix: Add `backend: TkAgg` to `~/.matplotlib/matplotlibrc` file
- Performance in Docker is way bette, especially JIT compilation is 4x faster
- Get inlcude paths: `echo | gcc -E -Wp,-v -`


### Path to use Bessel functions
C++17 functions cannpot be used. Boost functions also not per default. `Expression("boost::math::cyl_bessel_i(0,atan2(x[1], x[0]))", degree=2)` is allowed if one changes in file `/usr/local/lib/python3.6/dist-packages/dolfin/jit/jit.py`:

```
_math_header = """
// cmath functions
#include <boost/math/special_functions/bessel.hpp> // Added
%s
```

## Some Notes About Python

- Get current work directory:

```python
import os
cwd = os.getcwd()
print(cwd)
```

- Latex font for matplotlib

```python
# LaTeX text fonts:
# Use with raw strings: r"$\mathcal{O}(h^1)$"
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
```

- Get system path where modules are searched:

```
import sys
print(sys.path)
```

## Gitlab CI
- In `~/.gitlab-runner/config.toml`(for the runner):
  - change priviliges to ~true~
  - Use local images: `pull_policy = "if-not-present"`
- Run local: `gitlab-runner exec docker --docker-privileged build` or with `build` replaced by job name
  - maybe local vars have to be change to use local Docker images because `CI_REGISTRY`,... are not set
