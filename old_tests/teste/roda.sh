touch FD.in
nohup /opt/petsc/petsc/petsc_debug_0/bin/mpirun -n 2 ./MD3D_4.9_swarm -rtol 1.0E-6 -denok 1.0E-9 <FD.in >FD.out &

