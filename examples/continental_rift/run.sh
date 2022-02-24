/Users/victorsacek/Documents/petsc/arch-label-optimized/bin/mpirun -n 6  \
/Users/victorsacek/Documents/gits/md3d/MD3D_5.0_swarm \
-denok 1.0E-10 \
-particles_per_ele 100 \
-xi_min 1.0e-7 \
-ve 1 \
-te 1 \
-visc_const_per_element 1 \
-visc_harmonic_mean 0 \
-initial_dynamic_range 1 \
-seed 0,2 -strain_seed 0.0,1.0
