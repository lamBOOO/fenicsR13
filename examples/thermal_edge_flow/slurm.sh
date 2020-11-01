#!/usr/local_rwth/bin/zsh

#SBATCH --ntasks=100
#SBATCH --time=0-01:00:00
#SBATCH --mem-per-cpu=30000M
#SBATCH --job-name=fenicsR13_test
#SBATCH --output=output.%J.txt

cd $HOME/fenicsR13/examples/thermal_edge_flow
$MPIEXEC $FLAGS_MPI_BATCH singularity exec /rwthfs/rz/SW/UTIL.common/singularity/fenicsr13_latest.sif fenicsR13 input9.yml
