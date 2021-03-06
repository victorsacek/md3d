static char help[] = "S1\n";

/* Contributed by Dave May */

#include <petscksp.h>
#include <petscdmda.h>

#include <petsctime.h>

#include "petscsys.h"

#include "header.h"


PetscErrorCode create_thermal_3d(PetscInt mx,PetscInt my,PetscInt mz,PetscInt Px,PetscInt Py,PetscInt Pz);

PetscErrorCode build_thermal_3d();

PetscErrorCode solve_thermal_3d();

PetscErrorCode destroy_thermal_3d();

PetscErrorCode write_thermal_3d(int cont);

PetscErrorCode write_pressure(int cont);

PetscErrorCode write_geoq_3d(int cont);

PetscErrorCode create_veloc_3d(PetscInt mx,PetscInt my,PetscInt mz,PetscInt Px,PetscInt Py,PetscInt Pz);

PetscErrorCode createSwarm();
PetscErrorCode moveSwarm(PetscReal dt);
PetscErrorCode Swarm_add_remove();
PetscErrorCode SwarmViewGP(DM dms,const char prefix[]);



PetscErrorCode reader(int rank);

PetscErrorCode write_veloc_3d(int cont);
PetscErrorCode write_veloc_cond(int cont);

PetscErrorCode destroy_veloc_3d();

PetscErrorCode Calc_dt_calor();

PetscErrorCode write_tempo(int cont);

PetscErrorCode veloc_total();

PetscErrorCode calc_drho();





int main(int argc,char **args)
{
	PetscErrorCode ierr;
	char prefix[PETSC_MAX_PATH_LEN];
	PetscInt       Px,Py,Pz;
	
	ierr = PetscInitialize(&argc,&args,(char*)0,help);CHKERRQ(ierr);
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	seed = rank;
	
	reader(rank);
	
	PetscLogDouble Tempo1,Tempo2;
	
	PetscTime(&Tempo1);
	
	
	Px   = Py = Pz = PETSC_DECIDE;
	ierr = PetscOptionsGetInt(NULL,NULL,"-Px",&Px,NULL);CHKERRQ(ierr);
	Py   = Px; Pz = Px;
	ierr = PetscOptionsGetInt(NULL,NULL,"-Py",&Py,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-Pz",&Pz,NULL);CHKERRQ(ierr);

	if (n_interfaces>0){
		ierr = PetscCalloc1(n_interfaces, &seed_layer); CHKERRQ(ierr);
		seed_layer_size = n_interfaces;
		ierr = PetscOptionsGetIntArray(NULL,NULL,"-seed",seed_layer,&seed_layer_size,&seed_layer_set); CHKERRQ(ierr);

		ierr = PetscCalloc1(n_interfaces, &strain_seed_layer); CHKERRQ(ierr);
		strain_seed_layer_size = n_interfaces;
		ierr = PetscOptionsGetRealArray(NULL,NULL,"-strain_seed",strain_seed_layer,&strain_seed_layer_size,&strain_seed_layer_set); CHKERRQ(ierr);
		if (strain_seed_layer_set == PETSC_TRUE && seed_layer_set == PETSC_FALSE) {
			PetscPrintf(PETSC_COMM_WORLD,"Specify the seed layer with the flag -seed (required by -strain_seed)\n");
			exit(1);
		}
		if (strain_seed_layer_set == PETSC_TRUE && seed_layer_set == PETSC_TRUE && seed_layer_size != strain_seed_layer_size) {
			PetscPrintf(PETSC_COMM_WORLD,"Specify the same number of values in the list for flags -seed and -strain_seed\n");
			exit(1);
		}
		if (strain_seed_layer_set == PETSC_FALSE && seed_layer_set == PETSC_TRUE) {
			PetscPrintf(PETSC_COMM_WORLD,"Using default value '2.0' for -strain_seed (for all seed layers)\n");
			for (int k = 0; k < seed_layer_size; k++) {
				strain_seed_layer[k] = 2.0;
			}
		}
		PetscPrintf(PETSC_COMM_WORLD,"Number of seed layers: %d\n", seed_layer_size);
		for (int k = 0; k < seed_layer_size; k++) {
			PetscPrintf(PETSC_COMM_WORLD,"seed layer: %d - strain: %lf\n", seed_layer[k], strain_seed_layer[k]);
		}
		PetscPrintf(PETSC_COMM_WORLD,"\n");
	}

	ierr = PetscOptionsGetReal(NULL,NULL,"-random_initial_strain",&random_initial_strain,NULL);CHKERRQ(ierr);
	
	rtol = PETSC_DEFAULT;
	ierr = PetscOptionsGetReal(NULL,NULL,"-rtol",&rtol,NULL);CHKERRQ(ierr);
	
	temper_extern = 0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-te",&temper_extern,NULL);CHKERRQ(ierr);
	
	veloc_extern = 0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-ve",&veloc_extern,NULL);CHKERRQ(ierr);
	
	bcv_extern = 0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-bcve",&bcv_extern,NULL);CHKERRQ(ierr);

	ierr = PetscOptionsGetInt(NULL,NULL,"-initial_dynamic_range",&initial_dynamic_range,NULL);CHKERRQ(ierr);	
	
	denok_min = 1.0E-4;
	ierr = PetscOptionsGetReal(NULL,NULL,"-denok",&denok_min,NULL);CHKERRQ(ierr);

	Xi_min = 1.0E-14;
	ierr = PetscOptionsGetReal(NULL,NULL,"-xi_min",&Xi_min,NULL);CHKERRQ(ierr);
	
	print_visc = 0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-print_visc",&print_visc,NULL);CHKERRQ(ierr);

	ierr = PetscOptionsGetInt(NULL,NULL,"-particles_per_ele",&particles_per_ele,NULL);CHKERRQ(ierr);
	
	visc_const_per_element=0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-visc_const_per_element",&visc_const_per_element,NULL);CHKERRQ(ierr);
	
	visc_harmonic_mean=1;
	ierr = PetscOptionsGetInt(NULL,NULL,"-visc_harmonic_mean",&visc_harmonic_mean,NULL);CHKERRQ(ierr);
	
	dx_const = Lx/(Nx-1);
	dy_const = Ly/(Ny-1);
	dz_const = depth/(Nz-1);
	
	if (rank==0) printf("%lf %lf %lf\n",dx_const,dy_const,dz_const);
	
	
	ierr = create_thermal_3d(Nx-1,Ny-1,Nz-1,Px,Py,Pz);CHKERRQ(ierr);

	ierr = write_thermal_3d(-1);
	
	ierr = create_veloc_3d(Nx-1,Ny-1,Nz-1,Px,Py,Pz);CHKERRQ(ierr);

	if (geoq_on){
		PetscPrintf(PETSC_COMM_WORLD,"Swarm INICIO\n");
		ierr = createSwarm();CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"Swarm FIM\n");
	}
	
	calc_drho();
	
	//ierr = veloc_total(); CHKERRQ(ierr);
	if (visc_MAX>visc_MIN && initial_dynamic_range>0){
		double visc_contrast = PetscLog10Real(visc_MAX/visc_MIN);

		double visc_mean = PetscPowReal(10.0,PetscLog10Real(visc_MIN)+visc_contrast/2);

		int n_visc=0;

		visc_MIN_comp = visc_mean;
		visc_MAX_comp = visc_mean;

		PetscPrintf(PETSC_COMM_WORLD,"\n\nViscosity range: %.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

		ierr = veloc_total(); CHKERRQ(ierr);

		while ((visc_MIN_comp!=visc_MIN) && (visc_MAX_comp!=visc_MAX)){

			visc_MIN_comp = visc_mean*PetscPowReal(10.0,-n_visc*1.0);
			visc_MAX_comp = visc_mean*PetscPowReal(10.0,n_visc*1.0);

			if (visc_MIN_comp<visc_MIN) visc_MIN_comp=visc_MIN;
			if (visc_MAX_comp>visc_MAX) visc_MAX_comp=visc_MAX;

			PetscPrintf(PETSC_COMM_WORLD,"\n\nViscosity range: %.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

			ierr = veloc_total(); CHKERRQ(ierr);

			n_visc++;
		}
	}
	else {
		if (visc_MAX==visc_MIN) visc_MAX = visc_MIN*1.0001;  //avoiding the problem to the f2 in the denominator (Gerya...)
		visc_MIN_comp = visc_MIN;
		visc_MAX_comp = visc_MAX;
		ierr = veloc_total(); CHKERRQ(ierr);
	}
	
	
	PetscPrintf(PETSC_COMM_WORLD,"passou veloc_total\n");
	
	ierr = write_veloc_3d(tcont);
	ierr = write_veloc_cond(tcont);
	ierr = write_thermal_3d(tcont);
	ierr = write_pressure(tcont);
	ierr = write_geoq_3d(tcont);
	ierr = write_tempo(tcont);
	
	VecCopy(Veloc_fut,Veloc);
	
	PetscPrintf(PETSC_COMM_WORLD,"passou impressao\n");
	
	ierr = Calc_dt_calor();
	
	//float aux_le;
	
	
	
	for (tempo = dt_calor,tcont=1;tempo<=timeMAX && tcont<=stepMAX;tempo+=dt_calor, tcont++){
		
		
		
		ierr = build_thermal_3d();CHKERRQ(ierr);

		ierr = solve_thermal_3d();CHKERRQ(ierr);
		
		ierr = veloc_total(); CHKERRQ(ierr);
		
		if (geoq_on){
			for (PetscInt cont=0, max_cont=10;cont<max_cont; cont++){
				double fac = (1.0/max_cont)*(0.5+cont);
				//PetscPrintf(PETSC_COMM_WORLD,"%f %f\n",fac,(1.0-fac));
				//        x   ->   y
				VecCopy(Veloc,Veloc_weight);
				//			y			a		b		x
				VecAXPBY(Veloc_weight, fac, (1.0-fac),Veloc_fut); //y = a*x + b*y
				ierr = moveSwarm(dt_calor_sec/max_cont);
			}
			Swarm_add_remove();
			//exit(1);
		}
		
		if (tcont%print_step==0){
			ierr = write_thermal_3d(tcont);
			ierr = write_geoq_3d(tcont);
			ierr = write_veloc_3d(tcont);
			ierr = write_tempo(tcont);
			ierr = write_pressure(tcont);
			PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN-1,"step_%d",tcont);
			if (geoq_on){
				ierr = SwarmViewGP(dms,prefix);CHKERRQ(ierr);
			}
		}
		

		
		//if (rank==0) scanf("%f",&aux_le);
		
		//MPI_Barrier(PETSC_COMM_WORLD);
		
		
		ierr = Calc_dt_calor();
		
	}
	if (rank==0) printf("write\n");
	
	
	
	ierr = destroy_thermal_3d();CHKERRQ(ierr);
	
	destroy_veloc_3d();

	PetscTime(&Tempo2);
	
	//if (rank==0) printf("Tempo: %lf\n",Tempo2-Tempo1);
	
	ierr = PetscFinalize();
	return 0;
}


PetscErrorCode Calc_dt_calor(){
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	//////calc dt
	PetscInt ind_v_max,ind_v_min,ind_v_mod;
	PetscReal min_v,max_v,max_mod_v,dh_v_mod;
	
	VecMax(Veloc_fut,&ind_v_min,&max_v);
	VecMin(Veloc_fut,&ind_v_max,&min_v);
	//printf("max_v = %g\n",max_v);
	//printf("max_v = %g\n",min_v);
	max_mod_v = fabs(max_v);
	ind_v_mod = ind_v_max;
	if (max_mod_v<fabs(min_v)){
		max_mod_v = fabs(min_v);
		ind_v_mod = ind_v_min;
	}
	if (ind_v_mod%3==0) dh_v_mod = dx_const;
	if (ind_v_mod%3==1) dh_v_mod = dy_const;
	if (ind_v_mod%3==2) dh_v_mod = dz_const;
	if (rank==0) printf("dt = %g",(dh_v_mod/max_mod_v)/seg_per_ano);
	dt_calor = 0.2*(dh_v_mod/max_mod_v)/seg_per_ano;
	if (dt_calor>dt_MAX) dt_calor=dt_MAX;
	dt_calor_sec = dt_calor*seg_per_ano;
	////////fim calc dt
	
	PetscFunctionReturn(0);
	
}


PetscErrorCode write_tempo(int cont){

	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	PetscViewer viewer;
	
	char nome[100];
	
	sprintf(nome,"Tempo_%d.txt",cont);
	
	PetscReal aa[1];
	
	aa[0]=tempo;
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	PetscRealView(1,aa,viewer);
	PetscViewerDestroy(&viewer);
	
	if (rank==0) printf("Tempo: %lf\n",tempo);
	
	
	
	PetscFunctionReturn(0);
	
}

