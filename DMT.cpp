#include <petscksp.h>
#include <petscdmda.h>

#include <petsctime.h>

PetscErrorCode montaKeThermal_general(PetscReal *Ke, PetscReal *Me, PetscReal *Fe);

PetscErrorCode DMDAGetElementCorners(DM da,PetscInt *sx,PetscInt *sy,PetscInt *sz,PetscInt *mx,PetscInt *my,PetscInt *mz);

PetscErrorCode DMDAGetLocalElementSize(DM da,PetscInt *mxl,PetscInt *myl,PetscInt *mzl);

PetscErrorCode DMDAGetElementCorners(DM da,PetscInt *sx,PetscInt *sy,PetscInt *sz,PetscInt *mx,PetscInt *my,PetscInt *mz);

PetscErrorCode AssembleA_Thermal(Mat TA,DM thermal_da,PetscReal *TKe,PetscReal *TMe,PetscReal *TFe,
								 DM veloc_da, Vec Veloc_total);

PetscErrorCode AssembleF_Thermal(Vec F,DM thermal_da,PetscReal *TKe,PetscReal *TMe,PetscReal *TFe,
								 DM veloc_da, Vec Veloc_total);

PetscErrorCode Thermal_init(Vec F,DM thermal_da);

extern double Lx, Ly, depth;

extern PetscReal *NT;
extern PetscReal *NT_x;
extern PetscReal *NT_y;
extern PetscReal *NT_z;

extern long GaussQuad;

extern long T_NE;

extern PetscReal *TKe, *TCe, *TFe, *TCe_fut, *TMe, *Ttotal, *Ttotal_b;

extern PetscReal *T_vec_aux_ele;
extern PetscReal *T_vec_aux_ele_final;

extern PetscReal *v_vec_aux_ele;

extern long V_GT;
extern long V_GN;

extern Mat TA, TB;
extern Vec Tf, Temper;

extern Vec dRho;

extern KSP T_ksp;

extern DM da_Thermal;


extern Vec Veloc, Veloc_fut;


extern DM da_Veloc;


PetscErrorCode create_thermal_3d(PetscInt mx,PetscInt my,PetscInt mz,PetscInt Px,PetscInt Py,PetscInt Pz)
{

	PetscInt       dof,stencil_width;
	PetscErrorCode ierr;
	
	
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	PetscLogDouble Tempo1p,Tempo2p;
	
	PetscTime(&Tempo1p);

	
	PetscFunctionBeginUser;
	
	
	dof           = 1;
	stencil_width = 1;
	ierr          = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
								 mx+1,my+1,mz+1,Px,Py,Pz,dof,stencil_width,NULL,NULL,NULL,&da_Thermal);CHKERRQ(ierr);
	ierr = DMDASetFieldName(da_Thermal,0,"T");CHKERRQ(ierr);

	
	
	
	
	ierr = PetscCalloc1(GaussQuad*T_NE,&NT); CHKERRQ(ierr);
	ierr = PetscCalloc1(GaussQuad*T_NE,&NT_x); CHKERRQ(ierr);
	ierr = PetscCalloc1(GaussQuad*T_NE,&NT_y); CHKERRQ(ierr);
	ierr = PetscCalloc1(GaussQuad*T_NE,&NT_z); CHKERRQ(ierr);
	
	ierr = PetscCalloc1(T_NE*T_NE,&TKe); CHKERRQ(ierr);
	ierr = PetscCalloc1(T_NE*T_NE,&TCe); CHKERRQ(ierr);
	ierr = PetscCalloc1(T_NE*T_NE,&TCe_fut); CHKERRQ(ierr);
	ierr = PetscCalloc1(T_NE*T_NE,&TMe); CHKERRQ(ierr);
	ierr = PetscCalloc1(T_NE*T_NE,&Ttotal); CHKERRQ(ierr);
	ierr = PetscCalloc1(T_NE*T_NE,&Ttotal_b); CHKERRQ(ierr);
	
	ierr = PetscCalloc1(T_NE,&TFe); CHKERRQ(ierr);
	
	ierr = PetscCalloc1(T_NE,&T_vec_aux_ele); CHKERRQ(ierr);
	ierr = PetscCalloc1(T_NE,&T_vec_aux_ele_final); CHKERRQ(ierr);
	
	ierr = PetscCalloc1(V_GT,&v_vec_aux_ele); CHKERRQ(ierr);
	
	montaKeThermal_general(TKe,TMe,TFe);
	
	
	ierr = DMDASetUniformCoordinates(da_Thermal,0.0,Lx,0.0,Ly,-depth,0.0);CHKERRQ(ierr);
	
		
	/* Generate a matrix with the correct non-zero pattern of type AIJ. This will work in parallel and serial */
	ierr = DMSetMatType(da_Thermal,MATAIJ);CHKERRQ(ierr);
	ierr = DMCreateMatrix(da_Thermal,&TA);CHKERRQ(ierr);
	ierr = DMCreateMatrix(da_Thermal,&TB);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Thermal,&Temper);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Thermal,&Tf);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Thermal,&dRho);CHKERRQ(ierr);
	
	ierr = Thermal_init(Temper,da_Thermal);
	
	/*PetscViewer viewer;
	
	char nome[100];
	
	sprintf(nome,"Temper_0.txt");
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	VecView(Temper,viewer);
	PetscViewerDestroy(&viewer);*/
	
	
	
	ierr = KSPCreate(PETSC_COMM_WORLD,&T_ksp);CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(T_ksp,"thermal_"); /* stokes */ CHKERRQ(ierr);
	
	PetscTime(&Tempo2p);
	if (rank==0) printf("Thermal create : %lf\n",Tempo2p-Tempo1p);
	
	
	PetscFunctionReturn(0);
}


PetscErrorCode build_thermal_3d()
{

	PetscErrorCode ierr;
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	PetscLogDouble Tempo1p,Tempo2p;
	
	PetscTime(&Tempo1p);
	//if (rank==0) printf("TA,TB,Tf -> zero entries\n");
	ierr = MatZeroEntries(TA);CHKERRQ(ierr);
	//if (rank==0) printf("passou TA\n");
	ierr = MatZeroEntries(TB);CHKERRQ(ierr);
	ierr = VecZeroEntries(Tf);CHKERRQ(ierr);
	
	//if (rank==0) printf("build TA,Tf\n");
	ierr = AssembleA_Thermal(TA,da_Thermal,TKe,TMe,TFe,da_Veloc,Veloc_fut);CHKERRQ(ierr);
	//if (rank==0) printf("t\n");
	
	
	ierr = AssembleF_Thermal(Tf,da_Thermal,TKe,TMe,TFe,da_Veloc,Veloc);CHKERRQ(ierr);
	//if (rank==0) printf("t\n");

	
	PetscTime(&Tempo2p);
	if (rank==0) printf("Thermal build: %lf\n",Tempo2p-Tempo1p);
	
	PetscFunctionReturn(0);
		
}

PetscErrorCode solve_thermal_3d()
{
	PetscErrorCode ierr;
	PetscLogDouble Tempo1,Tempo2;
	
	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	PetscTime(&Tempo1);
	
	/* SOLVE */
	
	//if (rank==0) printf("k\n");
	ierr = KSPSetOperators(T_ksp,TA,TA);CHKERRQ(ierr);
	//if (rank==0) printf("k\n");
	ierr = KSPSetFromOptions(T_ksp);CHKERRQ(ierr);
	//if (rank==0) printf("k\n");
	
	ierr = KSPSetTolerances(T_ksp,1.0E-9,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	
	ierr = KSPSolve(T_ksp,Tf,Temper);CHKERRQ(ierr);

	PetscTime(&Tempo2);
	if (rank==0) printf("Thermal solve: %lf\n",Tempo2-Tempo1);
	
	PetscFunctionReturn(0);
	
}


PetscErrorCode destroy_thermal_3d()
{

	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	PetscLogDouble Tempo1,Tempo2;
	PetscTime(&Tempo1);
	
	PetscErrorCode ierr;
	ierr = KSPDestroy(&T_ksp);CHKERRQ(ierr);
	ierr = VecDestroy(&Temper);CHKERRQ(ierr);
	ierr = VecDestroy(&Tf);CHKERRQ(ierr);
	ierr = MatDestroy(&TA);CHKERRQ(ierr);
	ierr = MatDestroy(&TB);CHKERRQ(ierr);
	ierr = DMDestroy(&da_Thermal);CHKERRQ(ierr);
	
	PetscTime(&Tempo2);
	if (rank==0) printf("Thermal destroy: %lf\n",Tempo2-Tempo1);
	
	
	PetscFunctionReturn(0);	
}

PetscErrorCode write_thermal_3d(int cont)
{
	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	PetscLogDouble Tempo1,Tempo2;
	PetscTime(&Tempo1);
	
	PetscViewer viewer;
	
	char nome[100];
	
	sprintf(nome,"Temper_%d.txt",cont);
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	VecView(Temper,viewer);
	PetscViewerDestroy(&viewer);
	
	PetscTime(&Tempo2);
	if (rank==0) printf("Thermal write: %lf\n",Tempo2-Tempo1);
	
	PetscFunctionReturn(0);	
}

