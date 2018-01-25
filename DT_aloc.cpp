#include <petscksp.h>

extern Vec fh_DT;
extern Mat Kh_DT;

extern Vec h_DT;

extern PetscReal *hKe;

extern PetscInt *Indices_Kh_DT_I;
extern PetscInt *Indices_Kh_DT_J;

extern long n_quad_top;

extern PetscInt *cond_h_DT;

extern long nodes_top;

extern long Nx;
extern long Ny;

extern PetscInt *Quad_mesh;

extern PetscReal *Vfe_DT;
extern PetscInt *Indices_Vfe_DT;

extern KSP DT_ksp;

PetscErrorCode DT_aloc()
{
	PetscErrorCode ierr=0;
	
	printf("teste\n");
	
	nodes_top = Nx*Ny;
	n_quad_top = (Nx-1)*(Ny-1);

	ierr = VecCreate(PETSC_COMM_SELF,&h_DT);CHKERRQ(ierr);
	//ierr = VecCreateSeq(PETSC_COMM_SELF,nodes_top,&h_DT);CHKERRQ(ierr);
	
	
	ierr = PetscObjectSetName((PetscObject) h_DT, "h_DT");CHKERRQ(ierr);
	ierr = VecSetSizes(h_DT,PETSC_DECIDE,nodes_top);CHKERRQ(ierr);
	
	ierr = VecSetFromOptions(h_DT);CHKERRQ(ierr);
	ierr = VecDuplicate(h_DT,&fh_DT);CHKERRQ(ierr);
	
	
	ierr = PetscCalloc1(n_quad_top*4,&Quad_mesh); CHKERRQ(ierr);
	ierr = PetscCalloc1(4*4,&hKe); CHKERRQ(ierr);
	ierr = PetscCalloc1(4,&Vfe_DT); CHKERRQ(ierr);
	ierr = PetscCalloc1(4,&Indices_Vfe_DT); CHKERRQ(ierr);
	
	ierr = PetscCalloc1(4,&Indices_Kh_DT_I); CHKERRQ(ierr);
	ierr = PetscCalloc1(4,&Indices_Kh_DT_J); CHKERRQ(ierr);
	
	
	
	int cont,i,j;
	
	cont=0;
	for (i=0;i<Nx-1;i++){
		for (j=0;j<Ny-1;j++){
			Quad_mesh[cont*4+0] = (i)*Ny + j;
			Quad_mesh[cont*4+1] = (i)*Ny + j+1;
			Quad_mesh[cont*4+2] = (i+1)*Ny + j;
			Quad_mesh[cont*4+3] = (i+1)*Ny + j+1;
			cont++;
		}
	}
	printf("teste\n");
	
	ierr = MatCreate(PETSC_COMM_SELF,&Kh_DT);CHKERRQ(ierr);
	ierr = MatSetSizes(Kh_DT,PETSC_DECIDE,PETSC_DECIDE,nodes_top,nodes_top);CHKERRQ(ierr);
	ierr = MatSetFromOptions(Kh_DT);CHKERRQ(ierr);
	ierr = MatSetUp(Kh_DT);CHKERRQ(ierr);
	
	//ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nodes_top,nodes_top,20,0,&Kh_DT);
	
	KSPCreate(PETSC_COMM_SELF,&DT_ksp);
	
	return (ierr);
	
}