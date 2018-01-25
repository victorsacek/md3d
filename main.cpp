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

PetscErrorCode write_geoq_3d(int cont);

PetscErrorCode create_veloc_3d(PetscInt mx,PetscInt my,PetscInt mz,PetscInt Px,PetscInt Py,PetscInt Pz);

PetscErrorCode createSwarm();
PetscErrorCode moveSwarm(PetscReal dt);
PetscErrorCode Swarm_add_remove();
PetscErrorCode SwarmViewGP(DM dms,const char prefix[]);
PetscErrorCode Swarm2Mesh();

PetscErrorCode build_veloc_3d();

PetscErrorCode solve_veloc_3d();

PetscErrorCode reader(int rank);

PetscErrorCode write_veloc_3d(int cont);

PetscErrorCode destroy_veloc_3d();

PetscErrorCode Calc_dt_calor();

PetscErrorCode write_tempo(int cont);





int main(int argc,char **args)
{
	PetscErrorCode ierr;
	char prefix[PETSC_MAX_PATH_LEN];
	PetscInt       Px,Py,Pz;
	
	ierr = PetscInitialize(&argc,&args,(char*)0,help);CHKERRQ(ierr);
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	
	reader(rank);
	
	PetscLogDouble Tempo1,Tempo2;
	
	PetscTime(&Tempo1);
	
	
	Px   = Py = Pz = PETSC_DECIDE;
	ierr = PetscOptionsGetInt(NULL,NULL,"-Px",&Px,NULL);CHKERRQ(ierr);
	Py   = Px; Pz = Px;
	ierr = PetscOptionsGetInt(NULL,NULL,"-Py",&Py,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-Pz",&Pz,NULL);CHKERRQ(ierr);
	
	rtol = PETSC_DEFAULT;
	ierr = PetscOptionsGetReal(NULL,NULL,"-rtol",&rtol,NULL);CHKERRQ(ierr);
	
	temper_extern = 0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-te",&temper_extern,NULL);CHKERRQ(ierr);
	
	
	denok_min = 1.0E-4;
	ierr = PetscOptionsGetReal(NULL,NULL,"-denok",&denok_min,NULL);CHKERRQ(ierr);
	
	dx_const = Lx/(Nx-1);
	dy_const = Ly/(Ny-1);
	dz_const = depth/(Nz-1);
	
	if (rank==0) printf("%lf %lf %lf\n",dx_const,dy_const,dz_const);
	
	
	ierr = create_thermal_3d(Nx-1,Ny-1,Nz-1,Px,Py,Pz);CHKERRQ(ierr);
	
	ierr = create_veloc_3d(Nx-1,Ny-1,Nz-1,Px,Py,Pz);CHKERRQ(ierr);

	if (geoq_on){
		PetscPrintf(PETSC_COMM_WORLD,"Swarm INICIO\n");
		ierr = createSwarm();CHKERRQ(ierr);
		ierr = Swarm2Mesh();
		PetscPrintf(PETSC_COMM_WORLD,"Swarm FIM\n");
	}
	
	int tcont=0;
	
	
	ierr = build_veloc_3d();CHKERRQ(ierr);
	
	ierr = solve_veloc_3d();CHKERRQ(ierr);
	
	ierr = write_veloc_3d(tcont);
	ierr = write_thermal_3d(tcont);
	ierr = write_geoq_3d(tcont);
	ierr = write_tempo(tcont);
	
	VecCopy(Veloc_fut,Veloc);
	
	
	ierr = Calc_dt_calor();
	
	//float aux_le;
	
	
	
	for (tempo = dt_calor,tcont=1;tempo<=timeMAX && tcont<=stepMAX;tempo+=dt_calor, tcont++){
		
		
		
		ierr = build_thermal_3d();CHKERRQ(ierr);

		ierr = solve_thermal_3d();CHKERRQ(ierr);
		
		ierr = build_veloc_3d();CHKERRQ(ierr);
		
		ierr = solve_veloc_3d();CHKERRQ(ierr);
		
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
			ierr = Swarm2Mesh();
			//exit(1);
		}
		
		if (tcont%print_step==0){
			ierr = write_thermal_3d(tcont);
			ierr = write_geoq_3d(tcont);
			ierr = write_veloc_3d(tcont);
			ierr = write_tempo(tcont);
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

