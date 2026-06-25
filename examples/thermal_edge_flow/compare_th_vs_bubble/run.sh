for input in inputs/*.yml;
do
  mpirun -n 4 fenicsR13 "$input"
done
