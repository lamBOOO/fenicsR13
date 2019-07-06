# gasDynamicsFEM

Repository for Master thesis project regarding FEM simulations for non-equilibrium gas dynamics.

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
