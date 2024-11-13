# Thermal flow

## Execution

- Use `mesh8.h5` for good results
- Works with `mpirun -n 32` on `momentum @ ACoM`

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
