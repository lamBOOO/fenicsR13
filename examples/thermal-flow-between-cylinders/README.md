# Thermal flow

## Execution

- Use `mesh8.h5` for good results
- Adaptive meshes have better results
- Very sensitive to tau_p parameter with CIP, like thermal edge flow in Theisen2021
- Works with `mpirun -n 120` on `momentum @ ACoM`

## Post-processing

## Install Paraview with conda
```bash
conda create --name pv-test
conda activate pv-test
conda install conda-forge::paraview
# Alternatively, run the scripts using `pvpython`
# (located in e.g. `/Applications/ParaView-5.13.0.app/Contents/bin/pvpython` on macOS)
```

## VSCode linting
- Use the conda installation since `pvpython` doesn't seem to work properly.
