#include <petscksp.h>
#include <petscdmda.h>

#include <petsctime.h>


extern int bcv_top_normal;
extern int bcv_top_slip;

extern int bcv_bot_normal;
extern int bcv_bot_slip;

extern int bcv_left_normal;
extern int bcv_left_slip;

extern int bcv_right_normal;
extern int bcv_right_slip;


typedef struct {
	PetscScalar u;
	PetscScalar v;
	PetscScalar w;
	//PetscScalar p;
} Stokes;

PetscErrorCode AssembleA_Veloc(Mat A,Mat AG,DM veloc_da, DM temper_da);

PetscErrorCode AssembleF_Veloc(Vec F,DM veloc_da,DM drho_da, Vec FP);

PetscErrorCode montaKeVeloc_general(PetscReal *KeG, double dx_const, double dy_const, double dz_const);

PetscErrorCode montaKeVeloc_simplif(PetscReal *Ke,PetscReal *KeG,PetscReal *Temper_ele, PetscReal *geoq_ele);

PetscErrorCode montaCeVeloc(PetscReal *Ce);

PetscErrorCode montafeVeloc(PetscReal *fMe);

PetscErrorCode calc_drho();

PetscErrorCode write_veloc_3d(int cont);

double calc_visco_ponto(double T,double z,double geoq_ponto);

extern double r06;
extern double r8p9;
extern double r5p9;

extern long V_NE, V_GN, V_GT;

extern long GaussQuad;

extern Mat VA, VB, VG;
extern Vec Vf,Vf_P, Veloc, Veloc_fut, Veloc_weight;


extern Vec rk_vec2;

extern Vec rk_vec;
extern Vec sk_vec;
extern Vec gs_vec;
extern Vec uk_vec;

extern Vec zk_vec;
extern Vec zk_vec2;

extern Vec Veloc_Cond;

extern Vec Pressure;

extern DM da_Veloc;
extern DM da_Thermal;

extern KSP V_ksp;

extern double Lx, Ly, depth;

extern PetscReal *Ke_veloc;
extern PetscReal *Ke_veloc_final;
extern PetscReal *Ke_veloc_general;

extern PetscReal *VCe;
extern PetscReal *VfMe;

extern PetscReal *Vfe;


extern double dx_const;
extern double dy_const;
extern double dz_const;


extern Vec local_V;
extern Vec local_VC;
extern Vec local_FV;
extern Vec local_FP;
extern Vec local_P;

extern PetscReal rtol;

extern PetscReal denok_min;


extern Vec Precon;
extern Vec local_Precon;

extern double visc_aux_MAX;
extern double visc_aux_MIN;


PetscErrorCode create_veloc_3d(PetscInt mx,PetscInt my,PetscInt mz,PetscInt Px,PetscInt Py,PetscInt Pz)
{

	PetscInt       dof,stencil_width;

	PetscErrorCode ierr;
	
	
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	PetscLogDouble Tempo1p,Tempo2p;
	
	PetscTime(&Tempo1p);
	
	
	PetscFunctionBeginUser;
	
	
	dof           = 3; //modif
	stencil_width = 1;
	ierr          = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
								 mx+1,my+1,mz+1,Px,Py,Pz,dof,stencil_width,NULL,NULL,NULL,&da_Veloc);CHKERRQ(ierr);
	ierr = DMSetFromOptions(da_Veloc);CHKERRQ(ierr);
	ierr = DMSetUp(da_Veloc);CHKERRQ(ierr);
	
	ierr = DMDASetFieldName(da_Veloc,0,"V_x");CHKERRQ(ierr);
	ierr = DMDASetFieldName(da_Veloc,1,"V_y");CHKERRQ(ierr);
	ierr = DMDASetFieldName(da_Veloc,2,"V_z");CHKERRQ(ierr);
	//ierr = DMDASetFieldName(da_Veloc,3,"P");CHKERRQ(ierr);
	
	//printf("a\n");
	
	
	
	
	ierr = PetscCalloc1(V_GT*V_GT,&Ke_veloc); CHKERRQ(ierr);
	ierr = PetscCalloc1(V_GT*V_GT,&Ke_veloc_final); CHKERRQ(ierr);

	ierr = PetscCalloc1(V_GT*V_GT*GaussQuad,&Ke_veloc_general); CHKERRQ(ierr);
	//ierr = PetscCalloc1(V_GT,&Indices_Ke_veloc_I); CHKERRQ(ierr);
	//ierr = PetscCalloc1(V_GT,&Indices_Ke_veloc_J); CHKERRQ(ierr);
	
	ierr = PetscCalloc1(V_GT,&Vfe); CHKERRQ(ierr);
	//ierr = PetscCalloc1(V_GT,&Indices_Vfe); CHKERRQ(ierr);
	ierr = PetscCalloc1(V_GT*V_NE,&VfMe); CHKERRQ(ierr); // 24 x 8 no elemento hexa
	ierr = PetscCalloc1(V_GT,&VCe); CHKERRQ(ierr);
	
	
	ierr = DMDASetUniformCoordinates(da_Veloc,0.0,Lx,0.0,Ly,-depth,0.0);CHKERRQ(ierr);
	

	
	//printf("a\n");
	
	

	ierr = DMSetMatType(da_Veloc,MATAIJ);CHKERRQ(ierr);
	ierr = DMCreateMatrix(da_Veloc,&VA);CHKERRQ(ierr);
	ierr = DMCreateMatrix(da_Veloc,&VB);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&Veloc);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&Veloc_fut);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&Veloc_weight);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&Vf);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&Vf_P);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Veloc,&Pressure);CHKERRQ(ierr);
	
	ierr = DMCreateMatrix(da_Veloc,&VG);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Veloc,&rk_vec2);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Veloc,&rk_vec);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&sk_vec);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&gs_vec);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&uk_vec);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Veloc,&zk_vec);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&zk_vec2);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Veloc,&Precon);CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(da_Veloc,&Veloc_Cond);CHKERRQ(ierr);
	
	
	ierr = DMCreateLocalVector(da_Veloc,&local_V);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(da_Veloc,&local_VC);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(da_Veloc,&local_FV);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(da_Veloc,&local_FP);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(da_Veloc,&local_P);CHKERRQ(ierr);
	
	ierr = DMCreateLocalVector(da_Veloc,&local_Precon);CHKERRQ(ierr);
	
	
	Stokes					***ff;
	PetscInt               M,N,P;
	
	ierr = DMDAGetInfo(da_Veloc,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	
	
	ierr = VecZeroEntries(local_FV);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Veloc,local_FV,&ff);CHKERRQ(ierr);
	
	PetscInt       sx,sy,sz,mmx,mmy,mmz;
	PetscInt i,j,k;
	
	ierr = DMDAGetCorners(da_Veloc,&sx,&sy,&sz,&mmx,&mmy,&mmz);CHKERRQ(ierr);
	
	for (k=sz; k<sz+mmz; k++) {
		for (j=sy; j<sy+mmy; j++) {
			for (i=sx; i<sx+mmx; i++) {
				ff[k][j][i].u = 1.0;
				ff[k][j][i].v = 1.0;
				ff[k][j][i].w = 1.0;
				
				if (i==0   && bcv_left_normal==1) ff[k][j][i].u = 0.0;
				if (i==0   && bcv_left_slip==1) {
						ff[k][j][i].v = 0.0;
						ff[k][j][i].w = 0.0;
				}
				
				if (i==M-1 && bcv_right_normal==1)ff[k][j][i].u = 0.0;
				if (i==M-1   && bcv_right_slip==1) {
					ff[k][j][i].v = 0.0;
					ff[k][j][i].w = 0.0;
				}
				
				if (j==0 || j==N-1) ff[k][j][i].v = 0.0;
				
				if (k==0   && bcv_bot_normal==1) ff[k][j][i].w = 0.0;
				if (k==0   && bcv_bot_slip==1){
					ff[k][j][i].u = 0.0;
					ff[k][j][i].v = 0.0;
				}
				
				if (k==P-1 && bcv_top_normal==1) ff[k][j][i].w = 0.0;
				if (k==P-1 && bcv_top_slip==1){
					ff[k][j][i].u = 0.0;
					ff[k][j][i].v = 0.0;
				}
				
			}
		}
	}
	
	ierr = DMDAVecRestoreArray(da_Veloc,local_FV,&ff);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Veloc,local_FV,INSERT_VALUES,Veloc_Cond);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Veloc,local_FV,INSERT_VALUES,Veloc_Cond);CHKERRQ(ierr);
	
	int ind;
	PetscReal r;
	
	VecMax(Veloc_Cond,&ind,&r);
	
	if (rank==0) printf("VecMax Veloc_Cond: %f\n",r);
	
	VecSum(Veloc_Cond,&r);
	
	if (rank==0) printf("VecSum Veloc_Cond: %f\n",r);
	
	
	
	ierr = KSPCreate(PETSC_COMM_WORLD,&V_ksp);CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(V_ksp,"veloc_"); CHKERRQ(ierr);
	
	PetscTime(&Tempo2p);
	if (rank==0) printf("Velocity create : %lf\n",Tempo2p-Tempo1p);
	
	montaKeVeloc_general(Ke_veloc_general,dx_const,dy_const,dz_const);
	
	montaCeVeloc(VCe);
	montafeVeloc(VfMe);
	
	
	PetscFunctionReturn(0);
	 
}

PetscErrorCode build_veloc_3d()
{
	
	PetscErrorCode ierr;
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	PetscLogDouble Tempo1p,Tempo2p;
	
	PetscTime(&Tempo1p);
	if (rank==0) printf("VA,VB,Vf -> zero entries\n");
	ierr = MatZeroEntries(VA);CHKERRQ(ierr);
	if (rank==0) printf("passou VA\n");
	ierr = MatZeroEntries(VB);CHKERRQ(ierr);
	ierr = VecZeroEntries(Vf);CHKERRQ(ierr);
	ierr = VecZeroEntries(Vf_P);CHKERRQ(ierr);
	
	ierr = VecZeroEntries(Precon);CHKERRQ(ierr);
	
	if (rank==0) printf("build VA,Vf\n");
	ierr = AssembleA_Veloc(VA,VG,da_Veloc,da_Thermal);CHKERRQ(ierr);
	if (rank==0) printf("t\n");
	
	ierr = VecReciprocal(Precon);
	
	ierr = calc_drho();CHKERRQ(ierr);
	
	ierr = AssembleF_Veloc(Vf,da_Veloc,da_Thermal,Vf_P);CHKERRQ(ierr);
	if (rank==0) printf("t\n");
	
	
	PetscTime(&Tempo2p);
	if (rank==0) printf("Velocity build : %lf\n",Tempo2p-Tempo1p);
	
	PetscFunctionReturn(0);
	
}




PetscErrorCode solve_veloc_3d()
{
	PetscErrorCode ierr;
	PetscLogDouble Tempo1,Tempo2;
	
	int rank;
	
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	PetscTime(&Tempo1);
	
	
	VecCopy(Veloc_fut,Veloc);
	
	/* SOLVE */
	
	//if (rank==0) printf("k\n");
	ierr = KSPSetOperators(V_ksp,VA,VA);CHKERRQ(ierr);
	//if (rank==0) printf("k\n");
	ierr = KSPSetFromOptions(V_ksp);CHKERRQ(ierr);
	//if (rank==0) printf("k\n");
	ierr = KSPSetInitialGuessNonzero(V_ksp,PETSC_TRUE);
	
	ierr = KSPSetTolerances(V_ksp,rtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	
	
	////////
	
	PetscReal denok=0,betak,alphak;
	
	PetscInt maxk = 400,k;
	
	ierr = KSPSolve(V_ksp,Vf_P,Veloc_fut);CHKERRQ(ierr);
	
	//write_veloc_3d(101);
	
	ierr = MatMultTranspose(VG,Veloc_fut,rk_vec2);CHKERRQ(ierr);
	
	
	ierr = VecPointwiseMult(zk_vec2,Precon,rk_vec2);
	
	ierr = VecDot(rk_vec2,rk_vec2,&denok);CHKERRQ(ierr);
	
	if (rank==0) printf("denok = %lg\n",denok);
	
	for (k=1;k<maxk && denok>denok_min;k++){
		if (k==1) VecCopy(zk_vec2,sk_vec);
		else {
			VecCopy(zk_vec2,zk_vec);
			
			ierr = VecPointwiseMult(zk_vec2,Precon,rk_vec2);
			
			VecDot(zk_vec2,rk_vec2,&betak);
			VecDot(zk_vec,rk_vec,&denok);
			betak=betak/denok;
			VecAYPX(sk_vec,betak,zk_vec2);
		}
		ierr = MatMult(VG,sk_vec,gs_vec);
		
		VecPointwiseMult(gs_vec,gs_vec,Veloc_Cond);
		
		VecDot(gs_vec,gs_vec,&denok);
		
		
		
		KSPSolve(V_ksp,gs_vec,uk_vec);
		
		VecDot(zk_vec2,rk_vec2,&alphak);
		VecDot(gs_vec,uk_vec,&denok);
		
		alphak=alphak/denok;
		
		VecAXPY(Pressure,alphak,sk_vec);
		
		VecAXPY(Veloc_fut,-alphak,uk_vec);
		
		VecCopy(rk_vec2,rk_vec);
		
		ierr = MatMultTranspose(VG,uk_vec,rk_vec2);CHKERRQ(ierr);
		
		VecAYPX(rk_vec2,-alphak,rk_vec);
		
		VecDot(rk_vec2,rk_vec2,&denok);
		
		if (rank==0) printf("denok = %lg, k=%d\n",denok,k);
		
	}
	
	
	////////
	
	
	PetscTime(&Tempo2);
	if (rank==0) printf("Velocity solve: %lf\n",Tempo2-Tempo1);
	
	PetscFunctionReturn(0);
	
}


PetscErrorCode destroy_veloc_3d()
{
	
	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	PetscLogDouble Tempo1,Tempo2;
	PetscTime(&Tempo1);
	
	PetscErrorCode ierr;
	ierr = KSPDestroy(&V_ksp);CHKERRQ(ierr);
	ierr = VecDestroy(&Veloc_weight);CHKERRQ(ierr);
	ierr = VecDestroy(&Veloc_fut);CHKERRQ(ierr);
	ierr = VecDestroy(&Veloc);CHKERRQ(ierr);
	ierr = VecDestroy(&Vf);CHKERRQ(ierr);
	ierr = VecDestroy(&Vf_P);CHKERRQ(ierr);
	ierr = MatDestroy(&VA);CHKERRQ(ierr);
	ierr = MatDestroy(&VB);CHKERRQ(ierr);
	ierr = DMDestroy(&da_Veloc);CHKERRQ(ierr);
	
	PetscTime(&Tempo2);
	if (rank==0) printf("Velocity destroy: %lf\n",Tempo2-Tempo1);
	
	
	PetscFunctionReturn(0);
}

PetscErrorCode write_veloc_3d(int cont)
{
	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	PetscLogDouble Tempo1,Tempo2;
	PetscTime(&Tempo1);
	
	PetscViewer viewer;
	
	char nome[100];
	
	sprintf(nome,"Veloc_fut_%d.txt",cont);
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	VecView(Veloc_fut,viewer);
	PetscViewerDestroy(&viewer);
	
	PetscTime(&Tempo2);
	if (rank==0) printf("Velocity write: %lf\n",Tempo2-Tempo1);
	
	PetscFunctionReturn(0);
}










PetscErrorCode montaKeVeloc_general(PetscReal *KeG, double dx_const, double dy_const, double dz_const){
	
	long i,j,ii;
	long aux;
	
	
	
	double kx,ky,kz;
	
	double ex,ey,ez;
	
	long cont;
	
	double N_x[V_NE];
	double N_y[V_NE];
	double N_z[V_NE];
	
	double SN[6][V_GT];
	
	for (i=0;i<6;i++){
		for (j=0;j<V_GT;j++){
			SN[i][j]=0;
		}
	}
	
	
	long point;
	
	for (point=0;point<GaussQuad;point++){
		for (i=0;i<V_GT;i++){
			for (j=0;j<V_GT;j++){
				KeG[(i*V_GT+j)+point*V_GT*V_GT]=0.0;
			}
		}
	}
	
	double Hx,Hy,Hz,prodH;
	
	
	point=0;
	for (kz=-r06; kz<=r06; kz+=r06){
		if (kz==0) Hz=r8p9;
		else Hz=r5p9;
		for (ky=-r06; ky<=r06; ky+=r06){
			if (ky==0) Hy=r8p9;
			else Hy=r5p9;
			for (kx=-r06; kx<=r06; kx+=r06){
				if (kx==0) Hx=r8p9;
				else Hx=r5p9;
				
				prodH = Hx*Hy*Hz;
				cont=0;
				for (ez=-1.;ez<=1.;ez+=2.){
					for (ey=-1.;ey<=1.;ey+=2.){
						for (ex=-1.;ex<=1.;ex+=2.){
							//N[cont]=(1+ex*kx)*(1+ey*ky)*(1+ez*kz)/8.0;
							N_x[cont]=ex*(1+ey*ky)*(1+ez*kz)/4.0/dx_const;
							N_y[cont]=(1+ex*kx)*ey*(1+ez*kz)/4.0/dy_const;
							N_z[cont]=(1+ex*kx)*(1+ey*ky)*ez/4.0/dz_const;
							cont++;
						}
					}
				}
				
				
				for (j=0;j<V_NE;j++){
					aux = j*V_GN;
					SN[0][aux  ]=N_x[j];
					SN[1][aux+1]=N_y[j];
					SN[2][aux+2]=N_z[j];
					
					SN[3][aux  ]=N_y[j];SN[3][aux+1]=N_x[j];
					SN[4][aux+1]=N_z[j];SN[4][aux+2]=N_y[j];
					SN[5][aux  ]=N_z[j];					SN[5][aux+2]=N_x[j];
					
					//Nd[aux  ]=N_x[j];
					//Nd[aux+1]=N_y[j];
					//Nd[aux+2]=N_z[j];
				}
				
				
				for (i=0;i<V_GT;i++){
					for (j=0;j<V_GT;j++){
						for (ii=0;ii<3;ii++){
							KeG[(i*V_GT+j)+point*V_GT*V_GT]+=prodH*2*SN[ii][i]*SN[ii][j];
						}
						for (;ii<6;ii++){
							KeG[(i*V_GT+j)+point*V_GT*V_GT]+=prodH*SN[ii][i]*SN[ii][j];
						}
						
					}
				}
				
				point++;
				
			}
		}
	}
	
	
	PetscFunctionReturn(0);
	
}



PetscErrorCode montaKeVeloc_simplif(PetscReal *Ke,PetscReal *KeG,PetscReal *Temper_ele, PetscReal *geoq_ele){
	
	long i,j;
	
	double Visc_local,Temper_local,Geoq_local;
	
	
	double kx,ky,kz;
	
	double ex,ey,ez;
	
	long cont;
	
	//PetscReal Visc_ele[V_NE];
	
	//for (i=0;i<V_NE;i++) Visc_ele[i]=visco_const;
	
	//VecGetValues(visc_vec,V_NE,&Hexa_thermal[t*V_NE],Visc_ele);
	
	//for (i=0;i<V_NE;i++) printf("%g ",Visc_ele[i]);
	//printf("\n");
	
	double Hx,Hy,Hz,prodH;
	
	long point=0;
	
	for (i=0;i<V_GT*V_GT;i++) Ke[i]=0.0;
	
	for (kz=-r06; kz<=r06; kz+=r06){
		if (kz==0) Hz=r8p9;
		else Hz=r5p9;
		for (ky=-r06; ky<=r06; ky+=r06){
			if (ky==0) Hy=r8p9;
			else Hy=r5p9;
			for (kx=-r06; kx<=r06; kx+=r06){
				if (kx==0) Hx=r8p9;
				else Hx=r5p9;
				
				Temper_local = 0.0;
				Geoq_local = 0.0;
				
				prodH = Hx*Hy*Hz;
				cont=0;
				for (ez=-1.;ez<=1.;ez+=2.){
					for (ey=-1.;ey<=1.;ey+=2.){
						for (ex=-1.;ex<=1.;ex+=2.){
							Temper_local+=Temper_ele[cont]*(1+ex*kx)*(1+ey*ky)*(1+ez*kz)/8.0;
							Geoq_local+=geoq_ele[cont]*(1+ex*kx)*(1+ey*ky)*(1+ez*kz)/8.0;
							cont++;
						}
					}
				}
				
				Visc_local = calc_visco_ponto(Temper_local,0.0/*(z) mudar!!!!*/,Geoq_local);
				
				if (Visc_local<visc_aux_MIN) visc_aux_MIN=Visc_local;
				if (Visc_local>visc_aux_MAX) visc_aux_MAX=Visc_local;
				
				for (i=0;i<V_GT;i++){
					for (j=0;j<V_GT;j++){
						Ke[i*V_GT+j]+=KeG[(i*V_GT+j)+point]*Visc_local;
					}
				}
				
				point+=V_GT*V_GT;
			}
		}
	}
	
	//printf("Visc_min = %lg; Visc_max = %lg\n",Visc_min,Visc_max);
	
	PetscFunctionReturn(0);

}


PetscErrorCode montafeVeloc(PetscReal *fMe)
{
	long i,j;
	
	double kx,ky,kz;
	
	double ex,ey,ez;
	
	long cont;
	
	double N[V_NE];
	
	for (i=0;i<V_GT;i++){
		for (j=0;j<V_NE;j++){
			fMe[i*V_NE+j]=0.0;
		}
	}
	
	double Hx,Hy,Hz,prodH;
	for (kz=-r06; kz<=r06; kz+=r06){
		if (kz==0) Hz=r8p9;
		else Hz=r5p9;
		for (ky=-r06; ky<=r06; ky+=r06){
			if (ky==0) Hy=r8p9;
			else Hy=r5p9;
			for (kx=-r06; kx<=r06; kx+=r06){
				if (kx==0) Hx=r8p9;
				else Hx=r5p9;
				
				
				prodH = Hx*Hy*Hz;
				cont=0;
				for (ez=-1.;ez<=1.;ez+=2.){
					for (ey=-1.;ey<=1.;ey+=2.){
						for (ex=-1.;ex<=1.;ex+=2.){
							N[cont]=(1+ex*kx)*(1+ey*ky)*(1+ez*kz)/8.0;
							cont++;
						}
					}
				}
				
				
				for (i=0;i<V_NE;i++){
					for (j=0;j<V_NE;j++){
						fMe[(i*V_GN+2)*V_NE+j]+=prodH*N[i]*N[j];
					}
				}
				
				
				
			}
		}
	}
	
	for (i=0;i<V_GT*V_NE;i++){
		fMe[i]*=dx_const*dy_const*dz_const;
	}

	
	PetscFunctionReturn(0);
}


PetscErrorCode montaCeVeloc(PetscReal *Ce){
	
	long i;
	
	/*double dx,dy,dz;
	 
	 //
	 
	 dx = xyz_thermal[Hexa_thermal[t][2]][0]-xyz_thermal[Hexa_thermal[t][0]][0];
	 dy = xyz_thermal[Hexa_thermal[t][1]][1]-xyz_thermal[Hexa_thermal[t][0]][1];
	 dz = xyz_thermal[Hexa_thermal[t][4]][2]-xyz_thermal[Hexa_thermal[t][0]][2];
	 
	 //*/
	
	double kx,ky,kz;
	
	double ex,ey,ez;
	
	long cont;
	
	double N_x[V_NE];
	double N_y[V_NE];
	double N_z[V_NE];
	
	for (i=0;i<V_GT;i++){
		Ce[i]=0.0;
	}
	
	double Hx,Hy,Hz,prodH;
	
	for (kz=-r06; kz<=r06; kz+=r06){
		if (kz==0) Hz=r8p9;
		else Hz=r5p9;
		for (ky=-r06; ky<=r06; ky+=r06){
			if (ky==0) Hy=r8p9;
			else Hy=r5p9;
			for (kx=-r06; kx<=r06; kx+=r06){
				if (kx==0) Hx=r8p9;
				else Hx=r5p9;
				
				
				prodH = Hx*Hy*Hz;
				cont=0;
				for (ez=-1.;ez<=1.;ez+=2.){
					for (ey=-1.;ey<=1.;ey+=2.){
						for (ex=-1.;ex<=1.;ex+=2.){
							//N[cont]=(1+ex*kx)*(1+ey*ky)*(1+ez*kz)/8.0;
							N_x[cont]=ex*(1+ey*ky)*(1+ez*kz)/4.0/dx_const;
							N_y[cont]=(1+ex*kx)*ey*(1+ez*kz)/4.0/dy_const;
							N_z[cont]=(1+ex*kx)*(1+ey*ky)*ez/4.0/dz_const;
							cont++;
						}
					}
				}
				
				for (i=0;i<V_NE;i++){
					Ce[i*V_GN]+=prodH*N_x[i]; // 1/8 funcao de forma da pressao (constante no elemento)
					Ce[i*V_GN+1]+=prodH*N_y[i];
					Ce[i*V_GN+2]+=prodH*N_z[i];
				}
				
			}
		}
	}
	
	///adicionado
	for (i=0;i<V_GT;i++){
		Ce[i]*=-dx_const*dy_const*dz_const;
	}
	
	PetscFunctionReturn(0);

	
}