#include <petscksp.h>
#include <petscdmda.h>
#include <petscsys.h>
#include <petsctime.h>

PetscErrorCode DMDAGetElementCorners(DM da,PetscInt *sx,PetscInt *sy,PetscInt *sz,PetscInt *mx,PetscInt *my,PetscInt *mz);

PetscErrorCode montaKeVeloc_simplif(PetscReal *Ke,PetscReal *KeG, PetscReal visco_const);


typedef struct {
	PetscScalar u;
	PetscScalar v;
	PetscScalar w;
	//PetscScalar p;
} Stokes;


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

extern DM da_Veloc;

extern Vec dRho;

extern PetscInt Verif_VG;

extern double visco_r;

PetscErrorCode AssembleA_Veloc(Mat A,Mat AG,DM veloc_da){
	
	PetscErrorCode         ierr;
	
	PetscInt i,j,k;
	
	PetscInt               M,N,P;
	
	MatStencil ind1[1],ind[V_GT], indp[1];
	
	PetscScalar ONE[1]={1.0},u[V_GT*V_GT];
	
	//PetscScalar VCe_aux[V_GT];
	
	
	PetscInt       sx,sy,sz,mmx,mmy,mmz;
	
	ierr = DMDAGetCorners(veloc_da,&sx,&sy,&sz,&mmx,&mmy,&mmz);CHKERRQ(ierr);
	
	ierr = DMDAGetInfo(veloc_da,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	
	for (k=sz; k<sz+mmz; k++) {
		for (j=sy; j<sy+mmy; j++) {
			for (i=sx; i<sx+mmx; i++) {
				if (i==0 || i==M-1){
					ind1[0].i = i; ind1[0].j = j; ind1[0].k = k; ind1[0].c = 0;
					ierr = MatSetValuesStencil(A,1,ind1,1,ind1,ONE,ADD_VALUES);
				}
				if (j==0 || j==N-1){
					ind1[0].i = i; ind1[0].j = j; ind1[0].k = k; ind1[0].c = 1;
					ierr = MatSetValuesStencil(A,1,ind1,1,ind1,ONE,ADD_VALUES);
				}
				if (k==0 || k==P-1){
					ind1[0].i = i; ind1[0].j = j; ind1[0].k = k; ind1[0].c = 2;
					ierr = MatSetValuesStencil(A,1,ind1,1,ind1,ONE,ADD_VALUES);
				}
				/*if (i==M-1 || j==N-1 || k==P-1){
					ind1[0].i = i; ind1[0].j = j; ind1[0].k = k; ind1[0].c = 3;
					ierr = MatSetValuesStencil(A,1,ind1,1,ind1,ONE,ADD_VALUES);
				}*/
			}
		}
	}
	
	int n,g;
	
	PetscInt               sex,sey,sez,mx,my,mz;
	PetscInt               ei,ej,ek;
	
	ierr = DMDAGetElementCorners(veloc_da,&sex,&sey,&sez,&mx,&my,&mz);CHKERRQ(ierr);
	
	//printf("%d %d %d %d %d %d\n",sez,sez+mz-1,sey,sey+my-1,sex,sex+mx-1);
	
	PetscReal volume = dx_const*dy_const*dz_const;

	for (ek = sez; ek < sez+mz; ek++) {
		for (ej = sey; ej < sey+my; ej++) {
			for (ei = sex; ei < sex+mx; ei++) {
				
				montaKeVeloc_simplif(Ke_veloc,Ke_veloc_general,visco_r );
				
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
					if (ind[g].i==0 || ind[g].i==M-1){
						for (j=0;j<V_GT;j++) u[g*V_GT+j]=0.0;
					}
					else {
						for (j=0;j<V_GT;j++) u[g*V_GT+j]=Ke_veloc_final[g*V_GT+j];
					}
					
					g = 3*n+1;
					if (ind[g].j==0 || ind[g].j==N-1){
						for (j=0;j<V_GT;j++) u[g*V_GT+j]=0.0;
					}
					else {
						for (j=0;j<V_GT;j++) u[g*V_GT+j]=Ke_veloc_final[g*V_GT+j];
					}
					
					g = 3*n+2;
					if (ind[g].k==0 || ind[g].k==P-1){
						for (j=0;j<V_GT;j++) u[g*V_GT+j]=0.0;
					}
					else {
						for (j=0;j<V_GT;j++) u[g*V_GT+j]=Ke_veloc_final[g*V_GT+j];
					}
				}
				
				ierr = MatSetValuesStencil(A,V_GT,ind,V_GT,ind,u,ADD_VALUES);
				
				///
				
				//for (i=0;i<V_GT;i++) VCe_aux[i]=VCe[i];
				
				/*
				for (n=0;n<8;n++){
					g = 3*n;
					if (ind[g].i==0 || ind[g].i==M-1)	VCe_aux[g]=0.0;
					
					g = 3*n+1;
					if (ind[g].j==0 || ind[g].j==N-1)	VCe_aux[g]=0.0;
					
					g = 3*n+2;
					if (ind[g].k==0 || ind[g].k==P-1)	VCe_aux[g]=0.0;
				}*/
				//ierr = MatSetValuesStencil(A,V_GT,ind,1,indp,VCe_aux,ADD_VALUES);
				
				
				if (Verif_VG==0){
					indp[0].i=ei  ; indp[0].j=ej  ; indp[0].k=ek  ; indp[0].c=0;
					
					ierr = MatSetValuesStencil(AG,V_GT,ind,1,indp,VCe,ADD_VALUES);
				}
				
				
			}
		}
	}
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	
	if (Verif_VG==0){
		ierr = MatAssemblyBegin(AG,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(AG,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		
		Verif_VG=1;
	}
	
	/*char nome[100];
	PetscViewer viewer;
	sprintf(nome,"A_veloc.txt");
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	MatView(A,viewer);
	PetscViewerDestroy(&viewer);*/
	
	//printf("passou!!!\n");
	
	
	PetscFunctionReturn(0);
}

PetscErrorCode AssembleF_Veloc(Vec F,DM veloc_da,DM drho_da){
	
	Vec                    local_F,local_dRho;
	PetscScalar             ***rr;
	Stokes					***ff;
	
	
	PetscInt               M,N,P;
	PetscErrorCode         ierr;
	
	PetscScalar				Vfe[V_GT],dr[T_NE];
	
	PetscInt               sex,sey,sez,mx,my,mz;
	PetscInt               ei,ej,ek;
	
	
	PetscFunctionBeginUser;
	ierr = DMDAGetInfo(veloc_da,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	
	/* get acces to the vector */
	ierr = DMGetLocalVector(veloc_da,&local_F);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_F);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(veloc_da,local_F,&ff);CHKERRQ(ierr);
	
	
	ierr = DMGetLocalVector(drho_da,&local_dRho);CHKERRQ(ierr);
	
	ierr = VecZeroEntries(local_dRho);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(drho_da,dRho,INSERT_VALUES,local_dRho);
	ierr = DMGlobalToLocalEnd(  drho_da,dRho,INSERT_VALUES,local_dRho);
	
	ierr = DMDAVecGetArray(drho_da,local_dRho,&rr);CHKERRQ(ierr);
	
	PetscInt i,j,c,n,g;
	
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
					if (ind[g].i==0 || ind[g].i==M-1)	Vfe[g]=0.0;
					
					g = 3*n+1;
					if (ind[g].j==0 || ind[g].j==N-1)	Vfe[g]=0.0;
					
					g = 3*n+2;
					if (ind[g].k==0 || ind[g].k==P-1)	Vfe[g]=0.0;
				}
				
				
				for (c=0;c<V_GT;){
					ff[ind[c].k][ind[c].j][ind[c].i].u += Vfe[c]; c++;
					ff[ind[c].k][ind[c].j][ind[c].i].v += Vfe[c]; c++;
					ff[ind[c].k][ind[c].j][ind[c].i].w += Vfe[c]; c++;
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
	
	ierr = DMDAVecRestoreArray(veloc_da,local_F,&ff);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(veloc_da,local_F,ADD_VALUES,F);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(veloc_da,local_F,ADD_VALUES,F);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(veloc_da,&local_F);CHKERRQ(ierr);
	
	ierr = DMRestoreLocalVector(drho_da,&local_dRho);CHKERRQ(ierr);
	
	//printf("passou...\n");
	
	/*char nome[100];
	PetscViewer viewer;
	sprintf(nome,"F_veloc.txt");
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	VecView(F,viewer);
	PetscViewerDestroy(&viewer);*/
	
	PetscFunctionReturn(0);
}