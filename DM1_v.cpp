#include <petscksp.h>
#include <petscdmda.h>
#include <petscsys.h>
#include <petsctime.h>

PetscErrorCode DMDAGetElementCorners(DM da,PetscInt *sx,PetscInt *sy,PetscInt *sz,PetscInt *mx,PetscInt *my,PetscInt *mz);

PetscReal montaKeVeloc_simplif(PetscReal *Ke,PetscReal *KeG,PetscReal *Temper_ele, PetscReal *geoq_ele,
									PetscReal *e2_ele, PetscReal *VX_ele, PetscReal *VY_ele, PetscReal *VZ_ele,
									PetscReal z_base);

PetscErrorCode write_veloc_3d(int cont);
PetscErrorCode write_veloc_cond(int cont);

PetscErrorCode ascii2bin(char *s1, char *s2);

typedef struct {
	PetscScalar u;
	PetscScalar v;
	PetscScalar w;
} Stokes;

extern double Lx, Ly, depth;


extern PetscReal *VCe;

extern PetscReal *VfMe;

extern long V_GT;
extern long T_NE;

extern PetscReal *Ke_veloc;
extern PetscReal *Ke_veloc_final;
extern PetscReal *Ke_veloc_general;

extern double dx_const;
extern double dy_const;
extern double dz_const;

extern double depth;

extern DM da_Veloc;

extern Vec dRho;

extern Vec Temper;

extern PetscInt Verif_VG;

extern double visco_r;

extern Vec Veloc_Cond;

extern Vec Pressure;

extern Vec local_Temper;

extern Vec local_VC;

extern Vec local_FV;

extern Vec local_FP;

extern Vec local_P;

extern Vec local_dRho;

extern Vec geoq;
extern Vec local_geoq;

extern Vec geoq_strain;
extern Vec local_geoq_strain;


extern Vec Precon;
extern Vec local_Precon;


extern double visc_aux_MAX;
extern double visc_aux_MIN;

extern double e2_aux_MAX;
extern double e2_aux_MIN;

extern double seg_per_ano;

extern Vec Veloc;
extern Vec Veloc_fut;
extern Vec Veloc_0;

extern Vec local_V;

extern long Nx,Ny,Nz;

extern double depth;

extern PetscInt veloc_extern;

extern PetscReal print_visc;

extern int tcont;
extern long print_step;

extern Vec Adiag;

extern Mat VB;


PetscErrorCode AssembleA_Veloc(Mat A,Mat AG,DM veloc_da, DM temper_da){
	
	PetscErrorCode         ierr;
	
	PetscInt i,j,k;
	
	PetscInt               M,N,P;
	
	MatStencil indr[T_NE],ind1[1],ind[V_GT], indp[1];
	
	PetscScalar u[V_GT*V_GT],val_cond[1];
	
	Stokes					***pc;
	
	ierr = VecZeroEntries(local_Precon);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(veloc_da,local_Precon,&pc);CHKERRQ(ierr);
	
	/////
	Stokes					***VVC;
	//Vec						local_VC;
	
	ierr = VecZeroEntries(local_VC);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(veloc_da,Veloc_Cond,INSERT_VALUES,local_VC);
	ierr = DMGlobalToLocalEnd(  veloc_da,Veloc_Cond,INSERT_VALUES,local_VC);
	
	ierr = DMDAVecGetArray(veloc_da,local_VC,&VVC);CHKERRQ(ierr);
	//////
	
	PetscInt       sx,sy,sz,mmx,mmy,mmz;
	
	ierr = DMDAGetCorners(veloc_da,&sx,&sy,&sz,&mmx,&mmy,&mmz);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(veloc_da,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	
	for (k=sz; k<sz+mmz; k++) {
		for (j=sy; j<sy+mmy; j++) {
			for (i=sx; i<sx+mmx; i++) {
				
				ind1[0].i = i; ind1[0].j = j; ind1[0].k = k; ind1[0].c = 0;
				val_cond[0] = 1.0-VVC[k][j][i].u;
				ierr = MatSetValuesStencil(A,1,ind1,1,ind1,val_cond,ADD_VALUES);

				ind1[0].i = i; ind1[0].j = j; ind1[0].k = k; ind1[0].c = 1;
				val_cond[0] = 1.0-VVC[k][j][i].v;
				ierr = MatSetValuesStencil(A,1,ind1,1,ind1,val_cond,ADD_VALUES);
				
				ind1[0].i = i; ind1[0].j = j; ind1[0].k = k; ind1[0].c = 2;
				val_cond[0] = 1.0-VVC[k][j][i].w;
				ierr = MatSetValuesStencil(A,1,ind1,1,ind1,val_cond,ADD_VALUES);
				
				
			}
		}
	}
	
	int n,g;
	
	PetscInt               sex,sey,sez,mx,my,mz;
	PetscInt               ei,ej,ek;
	
	ierr = DMDAGetElementCorners(veloc_da,&sex,&sey,&sez,&mx,&my,&mz);CHKERRQ(ierr);
	
	//////////
	
	

	PetscScalar             ***tt;
	
	ierr = VecZeroEntries(local_Temper);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(temper_da,Temper,INSERT_VALUES,local_Temper);
	ierr = DMGlobalToLocalEnd(  temper_da,Temper,INSERT_VALUES,local_Temper);
	
	ierr = DMDAVecGetArray(temper_da,local_Temper,&tt);CHKERRQ(ierr);
	
	
	
	PetscScalar             ***qq;
	
	ierr = VecZeroEntries(local_geoq);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(temper_da,geoq,INSERT_VALUES,local_geoq);
	ierr = DMGlobalToLocalEnd(  temper_da,geoq,INSERT_VALUES,local_geoq);
	
	ierr = DMDAVecGetArray(temper_da,local_geoq,&qq);CHKERRQ(ierr);
	
	
	PetscScalar             ***qq_strain;
	
	ierr = VecZeroEntries(local_geoq_strain);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(temper_da,geoq_strain,INSERT_VALUES,local_geoq_strain);
	ierr = DMGlobalToLocalEnd(  temper_da,geoq_strain,INSERT_VALUES,local_geoq_strain);
	
	ierr = DMDAVecGetArray(temper_da,local_geoq_strain,&qq_strain);CHKERRQ(ierr);
	
	Stokes					***VV;
	
	ierr = VecZeroEntries(local_V);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(veloc_da,Veloc_fut,INSERT_VALUES,local_V);
	ierr = DMGlobalToLocalEnd(  veloc_da,Veloc_fut,INSERT_VALUES,local_V);
	
	ierr = DMDAVecGetArray(da_Veloc,local_V,&VV);CHKERRQ(ierr);
	
	
	
	PetscReal volume = dx_const*dy_const*dz_const;
	
	
	PetscReal temper_ele[T_NE],geoq_ele[T_NE],e2_ele[T_NE],VX_ele[T_NE],VY_ele[T_NE],VZ_ele[T_NE];
	
	visc_aux_MAX = 1.0E5;
	visc_aux_MIN = 1.0E50;
	
	e2_aux_MAX = 0.0;
	e2_aux_MIN = 1.0E50;
	
	PetscReal z_base;
	
	FILE *sai_visc;
	
	if (print_visc==1 && tcont%print_step==0) {
		
		char nome[200];
		
		PetscMPIInt rank;
		ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
		
		sprintf(nome,"visc_%d_%d.txt",tcont,rank);
		
		sai_visc = fopen(nome,"w");
		
	}
	
	PetscReal visc_meio;
	

	for (ek = sez; ek < sez+mz; ek++) {
		for (ej = sey; ej < sey+my; ej++) {
			for (ei = sex; ei < sex+mx; ei++) {
				
				
				indr[0].i=ei  ; indr[0].j=ej  ; indr[0].k=ek  ;
				indr[1].i=ei+1; indr[1].j=ej  ; indr[1].k=ek  ;
				indr[2].i=ei  ; indr[2].j=ej+1; indr[2].k=ek  ;
				indr[3].i=ei+1; indr[3].j=ej+1; indr[3].k=ek  ;
				indr[4].i=ei  ; indr[4].j=ej  ; indr[4].k=ek+1;
				indr[5].i=ei+1; indr[5].j=ej  ; indr[5].k=ek+1;
				indr[6].i=ei  ; indr[6].j=ej+1; indr[6].k=ek+1;
				indr[7].i=ei+1; indr[7].j=ej+1; indr[7].k=ek+1;
				
				for (i=0;i<T_NE;i++) temper_ele[i]=tt[indr[i].k][indr[i].j][indr[i].i];
				for (i=0;i<T_NE;i++) geoq_ele[i]=qq[indr[i].k][indr[i].j][indr[i].i];
				for (i=0;i<T_NE;i++) e2_ele[i]=qq_strain[indr[i].k][indr[i].j][indr[i].i];
				
				for (i=0;i<T_NE;i++) VX_ele[i]=VV[indr[i].k][indr[i].j][indr[i].i].u;
				for (i=0;i<T_NE;i++) VY_ele[i]=VV[indr[i].k][indr[i].j][indr[i].i].v;
				for (i=0;i<T_NE;i++) VZ_ele[i]=VV[indr[i].k][indr[i].j][indr[i].i].w;
				
				z_base = -depth + ek*dz_const;///!!!!
				
				
				visc_meio = montaKeVeloc_simplif(Ke_veloc,Ke_veloc_general,temper_ele,geoq_ele,e2_ele,VX_ele,VY_ele,VZ_ele,z_base);
				
				if (print_visc==1 && tcont%print_step==0) {
					fprintf(sai_visc,"%d %d %d %lg\n",ei,ej,ek,visc_meio);
				}
				
				for (i=0;i<V_GT*V_GT;i++) Ke_veloc_final[i]=Ke_veloc[i]*volume;
				
				n=0;
				ind[n].i=ei  ; ind[n].j=ej  ; ind[n].k=ek  ; ind[n].c=0; n++;
				ind[n].i=ei  ; ind[n].j=ej  ; ind[n].k=ek  ; ind[n].c=1; n++;
				ind[n].i=ei  ; ind[n].j=ej  ; ind[n].k=ek  ; ind[n].c=2; n++;
				
				ind[n].i=ei+1; ind[n].j=ej  ; ind[n].k=ek  ; ind[n].c=0; n++;
				ind[n].i=ei+1; ind[n].j=ej  ; ind[n].k=ek  ; ind[n].c=1; n++;
				ind[n].i=ei+1; ind[n].j=ej  ; ind[n].k=ek  ; ind[n].c=2; n++;
				
				ind[n].i=ei  ; ind[n].j=ej+1; ind[n].k=ek  ; ind[n].c=0; n++;
				ind[n].i=ei  ; ind[n].j=ej+1; ind[n].k=ek  ; ind[n].c=1; n++;
				ind[n].i=ei  ; ind[n].j=ej+1; ind[n].k=ek  ; ind[n].c=2; n++;
				
				ind[n].i=ei+1; ind[n].j=ej+1; ind[n].k=ek  ; ind[n].c=0; n++;
				ind[n].i=ei+1; ind[n].j=ej+1; ind[n].k=ek  ; ind[n].c=1; n++;
				ind[n].i=ei+1; ind[n].j=ej+1; ind[n].k=ek  ; ind[n].c=2; n++;
				
				ind[n].i=ei  ; ind[n].j=ej  ; ind[n].k=ek+1; ind[n].c=0; n++;
				ind[n].i=ei  ; ind[n].j=ej  ; ind[n].k=ek+1; ind[n].c=1; n++;
				ind[n].i=ei  ; ind[n].j=ej  ; ind[n].k=ek+1; ind[n].c=2; n++;
				
				ind[n].i=ei+1; ind[n].j=ej  ; ind[n].k=ek+1; ind[n].c=0; n++;
				ind[n].i=ei+1; ind[n].j=ej  ; ind[n].k=ek+1; ind[n].c=1; n++;
				ind[n].i=ei+1; ind[n].j=ej  ; ind[n].k=ek+1; ind[n].c=2; n++;
				
				ind[n].i=ei  ; ind[n].j=ej+1; ind[n].k=ek+1; ind[n].c=0; n++;
				ind[n].i=ei  ; ind[n].j=ej+1; ind[n].k=ek+1; ind[n].c=1; n++;
				ind[n].i=ei  ; ind[n].j=ej+1; ind[n].k=ek+1; ind[n].c=2; n++;
				
				ind[n].i=ei+1; ind[n].j=ej+1; ind[n].k=ek+1; ind[n].c=0; n++;
				ind[n].i=ei+1; ind[n].j=ej+1; ind[n].k=ek+1; ind[n].c=1; n++;
				ind[n].i=ei+1; ind[n].j=ej+1; ind[n].k=ek+1; ind[n].c=2; n++;
				
				
				for (n=0;n<8;n++){
					g = 3*n;
					if (VVC[indr[n].k][indr[n].j][indr[n].i].u==0){
						for (j=0;j<V_GT;j++) u[g*V_GT+j]=0.0;
					}
					else {
						for (j=0;j<V_GT;j++) u[g*V_GT+j]=Ke_veloc_final[g*V_GT+j];

					}
					
					g = 3*n+1;
					if (VVC[indr[n].k][indr[n].j][indr[n].i].v==0){
						for (j=0;j<V_GT;j++) u[g*V_GT+j]=0.0;
					}
					else {
						for (j=0;j<V_GT;j++) u[g*V_GT+j]=Ke_veloc_final[g*V_GT+j];
					}
					
					g = 3*n+2;
					if (VVC[indr[n].k][indr[n].j][indr[n].i].w==0){
						for (j=0;j<V_GT;j++) u[g*V_GT+j]=0.0;
					}
					else {
						for (j=0;j<V_GT;j++) u[g*V_GT+j]=Ke_veloc_final[g*V_GT+j];
					}
					
				}
				
				ierr = MatSetValuesStencil(A,V_GT,ind,V_GT,ind,u,ADD_VALUES);
				
				
				
				
				if (Verif_VG==0){
					indp[0].i=ei  ; indp[0].j=ej  ; indp[0].k=ek  ; indp[0].c=0;
					
					ierr = MatSetValuesStencil(AG,V_GT,ind,1,indp,VCe,ADD_VALUES);
				}
				
				
				//preconditioner construction
				for (n=0;n<V_GT;n++){
					pc[ek][ej][ei].u+=VCe[n]*VCe[n]/Ke_veloc_final[n*V_GT+n];
				}
				
				
			}
		}
	}
	
	if (print_visc==1 && tcont%print_step==0) {
		fclose(sai_visc);
	}
	
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	MatGetDiagonal(A,Adiag);
	VecReciprocal(Adiag);
	MatDiagonalSet(VB,Adiag,INSERT_VALUES);
	
	
	if (Verif_VG==0){
		ierr = MatAssemblyBegin(AG,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(AG,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		
		Verif_VG=1;
	}
	
	ierr = DMDAVecRestoreArray(da_Veloc,local_V,&VV);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(veloc_da,local_Precon,&pc);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(veloc_da,local_Precon,ADD_VALUES,Precon);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(veloc_da,local_Precon,ADD_VALUES,Precon);CHKERRQ(ierr);
	
	
	
	ierr = DMDAVecRestoreArray(veloc_da,local_VC,&VVC);
	ierr = DMDAVecRestoreArray(temper_da,local_Temper,&tt);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(temper_da,local_geoq,&qq);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(temper_da,local_geoq_strain,&qq_strain);CHKERRQ(ierr);
	
	printf("Visc_min = %lg, Visc_max = %lg\n",visc_aux_MIN,visc_aux_MAX);
	
	printf("e2_min = %lg, e2_max = %lg\n",e2_aux_MIN,e2_aux_MAX);
	
	PetscFunctionReturn(0);
}

PetscErrorCode AssembleF_Veloc(Vec F,DM veloc_da,DM drho_da,Vec FP){
	
	PetscScalar             ***rr;
	Stokes					***ff,***ffp,***VV;
	
	Stokes					***pp;
	
	
	PetscInt               M,N,P;
	PetscErrorCode         ierr;
	
	PetscScalar				Vfe[V_GT],Vfe_P[V_GT],dr[T_NE];
	
	PetscInt               sex,sey,sez,mx,my,mz;
	PetscInt               ei,ej,ek;
	
	
	PetscFunctionBeginUser;
	ierr = DMDAGetInfo(veloc_da,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	
	
	ierr = VecZeroEntries(local_V);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(veloc_da,Veloc,INSERT_VALUES,local_V);
	ierr = DMGlobalToLocalEnd(  veloc_da,Veloc,INSERT_VALUES,local_V);
	
	ierr = DMDAVecGetArray(veloc_da,local_V,&VV);CHKERRQ(ierr);
	
	
	
	/* get acces to the vector */

	ierr = VecZeroEntries(local_FV);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(veloc_da,local_FV,&ff);CHKERRQ(ierr);
	

	ierr = VecZeroEntries(local_FP);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(veloc_da,local_FP,&ffp);CHKERRQ(ierr);
	
	/////
	

	
	ierr = VecZeroEntries(local_dRho);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(drho_da,dRho,INSERT_VALUES,local_dRho);
	ierr = DMGlobalToLocalEnd(  drho_da,dRho,INSERT_VALUES,local_dRho);
	
	ierr = DMDAVecGetArray(drho_da,local_dRho,&rr);CHKERRQ(ierr);
	
	
	/////
	
	ierr = VecZeroEntries(local_P);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(veloc_da,Pressure,INSERT_VALUES,local_P);
	ierr = DMGlobalToLocalEnd(  veloc_da,Pressure,INSERT_VALUES,local_P);
	
	ierr = DMDAVecGetArray(veloc_da,local_P,&pp);CHKERRQ(ierr);
	
	/////
	Stokes					***VVC;
	//Vec						local_VC;
	
	ierr = VecZeroEntries(local_VC);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(veloc_da,Veloc_Cond,INSERT_VALUES,local_VC);
	ierr = DMGlobalToLocalEnd(  veloc_da,Veloc_Cond,INSERT_VALUES,local_VC);
	
	ierr = DMDAVecGetArray(veloc_da,local_VC,&VVC);CHKERRQ(ierr);
	//////
	
	

	
	
	
	PetscInt i,j,k,c,n,g;
	
	MatStencil indr[T_NE],ind[V_GT];
	
	ierr = DMDAGetElementCorners(drho_da,&sex,&sey,&sez,&mx,&my,&mz);CHKERRQ(ierr);
	
	PetscInt               sex1,sey1,sez1,mx1,my1,mz1;

	ierr = DMDAGetElementCorners(veloc_da,&sex1,&sey1,&sez1,&mx1,&my1,&mz1);CHKERRQ(ierr);
	
	if ((sex-sex1!=0)||(sey-sey1!=0)||(sez-sez1!=0)||(mx-mx1!=0)||(my-my1!=0)||(mz-mz1!=0)){
		printf("%d %d %d %d %d %d\n",sex-sex1,sey-sey1,sez-sez1,mx-mx1,my-my1,mz-mz1);
		SETERRQ1(PETSC_COMM_WORLD,1,"Wrong partition (temper,velocity)\n",1);
	}
	
	for (ek = sez; ek < sez+mz; ek++) {
		for (ej = sey; ej < sey+my; ej++) {
			for (ei = sex; ei < sex+mx; ei++) {
				//montaKeThermal_simplif(TCe, TKe, 0, 0);//modificar
				
				//for (c=0;c<T_NE*T_NE;c++) Ttotal_b[c] = TMe[c] - comp_alpha_thermal*dt_calor_sec*TCe[c];
				
				
				indr[0].i=ei  ; indr[0].j=ej  ; indr[0].k=ek  ;
				indr[1].i=ei+1; indr[1].j=ej  ; indr[1].k=ek  ;
				indr[2].i=ei  ; indr[2].j=ej+1; indr[2].k=ek  ;
				indr[3].i=ei+1; indr[3].j=ej+1; indr[3].k=ek  ;
				indr[4].i=ei  ; indr[4].j=ej  ; indr[4].k=ek+1;
				indr[5].i=ei+1; indr[5].j=ej  ; indr[5].k=ek+1;
				indr[6].i=ei  ; indr[6].j=ej+1; indr[6].k=ek+1;
				indr[7].i=ei+1; indr[7].j=ej+1; indr[7].k=ek+1;
				
				
				
				for (i=0;i<T_NE;i++) dr[i] = rr[indr[i].k][indr[i].j][indr[i].i];
				
				
				for (i=0;i<V_GT;i++){
					Vfe[i] = 0.0;
					for (j=0;j<T_NE;j++){
						Vfe[i]+=VfMe[i*T_NE+j]*dr[j];
					}
					Vfe_P[i] = Vfe[i] - VCe[i]*pp[ek][ej][ei].u;
				}
				
				n=0;
				ind[n].i=ei  ; ind[n].j=ej  ; ind[n].k=ek  ; ind[n].c=0; n++;
				ind[n].i=ei  ; ind[n].j=ej  ; ind[n].k=ek  ; ind[n].c=1; n++;
				ind[n].i=ei  ; ind[n].j=ej  ; ind[n].k=ek  ; ind[n].c=2; n++;
				
				ind[n].i=ei+1; ind[n].j=ej  ; ind[n].k=ek  ; ind[n].c=0; n++;
				ind[n].i=ei+1; ind[n].j=ej  ; ind[n].k=ek  ; ind[n].c=1; n++;
				ind[n].i=ei+1; ind[n].j=ej  ; ind[n].k=ek  ; ind[n].c=2; n++;
				
				ind[n].i=ei  ; ind[n].j=ej+1; ind[n].k=ek  ; ind[n].c=0; n++;
				ind[n].i=ei  ; ind[n].j=ej+1; ind[n].k=ek  ; ind[n].c=1; n++;
				ind[n].i=ei  ; ind[n].j=ej+1; ind[n].k=ek  ; ind[n].c=2; n++;
				
				ind[n].i=ei+1; ind[n].j=ej+1; ind[n].k=ek  ; ind[n].c=0; n++;
				ind[n].i=ei+1; ind[n].j=ej+1; ind[n].k=ek  ; ind[n].c=1; n++;
				ind[n].i=ei+1; ind[n].j=ej+1; ind[n].k=ek  ; ind[n].c=2; n++;
				
				ind[n].i=ei  ; ind[n].j=ej  ; ind[n].k=ek+1; ind[n].c=0; n++;
				ind[n].i=ei  ; ind[n].j=ej  ; ind[n].k=ek+1; ind[n].c=1; n++;
				ind[n].i=ei  ; ind[n].j=ej  ; ind[n].k=ek+1; ind[n].c=2; n++;
				
				ind[n].i=ei+1; ind[n].j=ej  ; ind[n].k=ek+1; ind[n].c=0; n++;
				ind[n].i=ei+1; ind[n].j=ej  ; ind[n].k=ek+1; ind[n].c=1; n++;
				ind[n].i=ei+1; ind[n].j=ej  ; ind[n].k=ek+1; ind[n].c=2; n++;
				
				ind[n].i=ei  ; ind[n].j=ej+1; ind[n].k=ek+1; ind[n].c=0; n++;
				ind[n].i=ei  ; ind[n].j=ej+1; ind[n].k=ek+1; ind[n].c=1; n++;
				ind[n].i=ei  ; ind[n].j=ej+1; ind[n].k=ek+1; ind[n].c=2; n++;
				
				ind[n].i=ei+1; ind[n].j=ej+1; ind[n].k=ek+1; ind[n].c=0; n++;
				ind[n].i=ei+1; ind[n].j=ej+1; ind[n].k=ek+1; ind[n].c=1; n++;
				ind[n].i=ei+1; ind[n].j=ej+1; ind[n].k=ek+1; ind[n].c=2; n++;
				
				
				
				for (n=0;n<8;n++){
					g = 3*n;
					if (VVC[indr[n].k][indr[n].j][indr[n].i].u==0){
						Vfe[g]=0.0;
						Vfe_P[g]=0.0;
					}
					
					g = 3*n+1;
					if (VVC[indr[n].k][indr[n].j][indr[n].i].v==0){
						Vfe[g]=0.0;
						Vfe_P[g]=0.0;
					}
					
					g = 3*n+2;
					if (VVC[indr[n].k][indr[n].j][indr[n].i].w==0){
						Vfe[g]=0.0;
						Vfe_P[g]=0.0;
					}
				}
				
				
				for (c=0;c<V_GT;){
					ff[ind[c].k][ind[c].j][ind[c].i].u += Vfe[c]; c++;
					ff[ind[c].k][ind[c].j][ind[c].i].v += Vfe[c]; c++;
					ff[ind[c].k][ind[c].j][ind[c].i].w += Vfe[c]; c++;
				}
				
				for (c=0;c<V_GT;){
					ffp[ind[c].k][ind[c].j][ind[c].i].u += Vfe_P[c]; c++;
					ffp[ind[c].k][ind[c].j][ind[c].i].v += Vfe_P[c]; c++;
					ffp[ind[c].k][ind[c].j][ind[c].i].w += Vfe_P[c]; c++;
				}
				
				
				/*
				for (c=0;c<T_NE;c++){
					T_vec_aux_ele_final[c]=dt_calor_sec*TFe[c];
					for (j=0;j<T_NE;j++){
						T_vec_aux_ele_final[c]+=tt[ind[j].k][ind[j].j][ind[j].i]*Ttotal_b[c*T_NE+j];
					}
				}
				
				for (i=0;i<8;i++){
					if (ind[i].k==0 || ind[i].k==P-1){
						T_vec_aux_ele_final[i]=0.0;
					}
				}
				
				for (c=0;c<T_NE;c++){
					ff[ind[c].k][ind[c].j][ind[c].i] += T_vec_aux_ele_final[c];
				}
				 */
			}
		}
	}
	
	
	PetscInt       sx,sy,sz,mmx,mmy,mmz;
	
	ierr = DMDAGetCorners(veloc_da,&sx,&sy,&sz,&mmx,&mmy,&mmz);CHKERRQ(ierr);
	
	for (k=sz; k<sz+mmz; k++) {
		for (j=sy; j<sy+mmy; j++) {
			for (i=sx; i<sx+mmx; i++) {
				if (VVC[k][j][i].u==0){
					ffp[k][j][i].u=VV[k][j][i].u;
				}
				if (VVC[k][j][i].v==0){
					ffp[k][j][i].v=VV[k][j][i].v;
				}
				if (VVC[k][j][i].w==0){
					ffp[k][j][i].w=VV[k][j][i].w;
				}
			}
		}
	}
	
	
	ierr = DMDAVecRestoreArray(veloc_da,local_FV,&ff);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(veloc_da,local_FV,ADD_VALUES,F);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(veloc_da,local_FV,ADD_VALUES,F);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(veloc_da,local_FP,&ffp);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(veloc_da,local_FP,ADD_VALUES,FP);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(veloc_da,local_FP,ADD_VALUES,FP);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(drho_da,local_dRho,&rr);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(veloc_da,local_P,&pp);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(veloc_da,local_VC,&VVC);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(veloc_da,local_V,&VV);CHKERRQ(ierr);
	
	//printf("passou...\n");
	
	/*char nome[100];
	PetscViewer viewer;
	sprintf(nome,"F_veloc.txt");
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	VecView(F,viewer);
	PetscViewerDestroy(&viewer);*/
	
	PetscFunctionReturn(0);
}

PetscErrorCode Init_Veloc(){
	
	
	PetscErrorCode         ierr;
	
	PetscMPIInt rank;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	
	
	
	PetscInt ix[1];
	PetscScalar y[1];
	
	PetscInt low,high;
	
	PetscInt t;

	
	if (veloc_extern==1){
		char s1[100],s2[100];
		
		sprintf(s1,"veloc_0_3D.txt");
		sprintf(s2,"veloc_init.bin");
		
		if (rank==0){
			ierr = ascii2bin(s1,s2); CHKERRQ(ierr);
		}
		MPI_Barrier(PETSC_COMM_WORLD);
		
		
		PetscInt size0;
		PetscViewer    viewer;
		
		VecGetSize(Veloc,&size0);
		
		
		Vec Fprov;
		
		//PetscPrintf(PETSC_COMM_WORLD,"size = %d\n",size);
		
		//PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Temper_init.bin",FILE_MODE_READ,&viewer);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD,s2,FILE_MODE_READ,&viewer);
		VecCreate(PETSC_COMM_WORLD,&Fprov);
		VecLoad(Fprov,viewer);
		PetscViewerDestroy(&viewer);
		
		
		VecGetOwnershipRange(Fprov,&low,&high);
		
		printf("%d %d\n",low,high);
		
		
		Vec FN;
		
		DMDACreateNaturalVector(da_Veloc,&FN);
		
		
		for (t=low;t<high;t++){
			ix[0] = t;
			VecGetValues(Fprov,1,ix,y);
			VecSetValue(FN,t,y[0], INSERT_VALUES);
		}
		
		VecAssemblyBegin(FN);
		VecAssemblyEnd(FN);
		
		DMDANaturalToGlobalBegin(da_Veloc,FN,INSERT_VALUES,Veloc);
		DMDANaturalToGlobalEnd(da_Veloc,FN,INSERT_VALUES,Veloc);
		
		
		//VecView(F,PETSC_VIEWER_STDOUT_WORLD);
		
		PetscBarrier(NULL);
		
		VecAssemblyBegin(Veloc);
		VecAssemblyEnd(Veloc);

	}
	else {
		
		Stokes					***VV;
		
		ierr = VecZeroEntries(local_V);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da_Veloc,local_V,&VV);CHKERRQ(ierr);
		
		PetscInt       sx,sy,sz,mmx,mmy,mmz;
		
		ierr = DMDAGetCorners(da_Veloc,&sx,&sy,&sz,&mmx,&mmy,&mmz);CHKERRQ(ierr);
		
		PetscInt               M,N,P;
		
		ierr = DMDAGetInfo(da_Veloc,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
		
		ierr = DMDAVecRestoreArray(da_Veloc,local_V,&VV);CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(da_Veloc,local_V,INSERT_VALUES,Veloc);CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd(da_Veloc,local_V,INSERT_VALUES,Veloc);CHKERRQ(ierr);
	}
	
	
	
	VecCopy(Veloc,Veloc_fut);
	VecCopy(Veloc,Veloc_0);
	
	write_veloc_3d(-1);
	write_veloc_cond(-1);
	
	PetscFunctionReturn(0);
	
}
