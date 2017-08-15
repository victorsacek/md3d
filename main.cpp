static char help[] = "S1\n";

/* Contributed by Dave May */

#include <petscksp.h>
#include <petscdmda.h>

#include <petsctime.h>

#include "header.h"


PetscErrorCode create_thermal_3d(PetscInt mx,PetscInt my,PetscInt mz,PetscInt Px,PetscInt Py,PetscInt Pz);

PetscErrorCode build_thermal_3d();

PetscErrorCode solve_thermal_3d();

PetscErrorCode destroy_thermal_3d();

PetscErrorCode write_thermal_3d(int cont);

PetscErrorCode create_veloc_3d(PetscInt mx,PetscInt my,PetscInt mz,PetscInt Px,PetscInt Py,PetscInt Pz);

PetscErrorCode build_veloc_3d();

PetscErrorCode solve_veloc_3d();

PetscErrorCode reader(int rank);

PetscErrorCode write_veloc_3d(int cont);

PetscErrorCode destroy_veloc_3d();

PetscErrorCode Calc_dt_calor();



int main(int argc,char **args)
{
	PetscErrorCode ierr;
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
	
	dx_const = Lx/(Nx-1);
	dy_const = Ly/(Ny-1);
	dz_const = depth/(Nz-1);
	
	if (rank==0) printf("%lf %lf %lf\n",dx_const,dy_const,dz_const);
	
	
	ierr = create_thermal_3d(Nx-1,Ny-1,Nz-1,Px,Py,Pz);CHKERRQ(ierr);
	
	ierr = create_veloc_3d(Nx-1,Ny-1,Nz-1,Px,Py,Pz);CHKERRQ(ierr);
	
	
	int tcont=0;
	
	
	ierr = build_veloc_3d();CHKERRQ(ierr);
	
	ierr = solve_veloc_3d();CHKERRQ(ierr);
	
	ierr = write_veloc_3d(tcont);
	ierr = write_thermal_3d(tcont);
	
	
	ierr = Calc_dt_calor();
	
	for (tempo = dt_calor;tempo<timeMAX && tcont<stepMAX;tempo+=dt_calor, tcont++){
		
		
		
		ierr = build_thermal_3d();CHKERRQ(ierr);

		ierr = solve_thermal_3d();CHKERRQ(ierr);
		
		
		ierr = build_veloc_3d();CHKERRQ(ierr);
		
		ierr = solve_veloc_3d();CHKERRQ(ierr);
		
		if (tcont%print_step==0){
			ierr = write_thermal_3d(tcont);
			ierr = write_veloc_3d(tcont);
		}
		
		
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