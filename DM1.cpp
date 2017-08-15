#include <petscksp.h>
#include <petscdmda.h>

#include <petsctime.h>

PetscErrorCode montaKeThermal_simplif(double *Ke_local,double *Ke);

extern PetscReal *TCe;
extern PetscReal *TCe_fut;

extern long T_NE;
extern PetscReal *Ttotal;
extern PetscReal *Ttotal_b;

extern double alpha_thermal;
extern double comp_alpha_thermal;
extern double dt_calor_sec;

extern PetscReal *T_vec_aux_ele;
extern PetscReal *T_vec_aux_ele_final;

extern PetscReal *v_vec_aux_ele;

extern Vec Temper;

extern Vec Temper_Cond;

extern double Delta_T;

extern int T_initial_cond;

typedef struct {
	PetscScalar u;
	PetscScalar v;
	PetscScalar w;
	//PetscScalar p;
} Stokes;


PetscErrorCode DMDAGetLocalElementSize(DM da,PetscInt *mxl,PetscInt *myl,PetscInt *mzl)
{
	PetscInt       m,n,p,M,N,P;
	PetscInt       sx,sy,sz;
	PetscErrorCode ierr;
	
	PetscFunctionBeginUser;
	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&sx,&sy,&sz,&m,&n,&p);CHKERRQ(ierr);
	
	if (mxl) {
		*mxl = m;
		if ((sx+m) == M) *mxl = m-1;  /* last proc */
	}
	if (myl) {
		*myl = n;
		if ((sy+n) == N) *myl = n-1;  /* last proc */
	}
	if (mzl) {
		*mzl = p;
		if ((sz+p) == P) *mzl = p-1;  /* last proc */
	}
	PetscFunctionReturn(0);
}


PetscErrorCode DMDAGetElementCorners(DM da,PetscInt *sx,PetscInt *sy,PetscInt *sz,PetscInt *mx,PetscInt *my,PetscInt *mz)
{
	PetscInt       si,sj,sk;
	PetscErrorCode ierr;
	
	PetscFunctionBeginUser;
	ierr = DMDAGetGhostCorners(da,&si,&sj,&sk,0,0,0);CHKERRQ(ierr);
	
	if (sx) {
		*sx = si;
		if (si != 0) *sx = si+1;
	}
	if (sy) {
		*sy = sj;
		if (sj != 0) *sy = sj+1;
	}
	if (sz) {
		*sz = sk;
		if (sk != 0) *sz = sk+1;
	}
	ierr = DMDAGetLocalElementSize(da,mx,my,mz);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


PetscErrorCode AssembleA_Thermal(Mat A,DM thermal_da,PetscReal *TKe,PetscReal *TMe,PetscReal *TFe,
								 DM veloc_da, Vec Veloc_total)
{
	//DM                     cda;
	
	PetscErrorCode ierr;
	
	PetscInt               sex,sey,sez,mx,my,mz;
	PetscInt               ei,ej,ek;
	
	
	Stokes					***VV;
	Vec						local_V;
	
	ierr = DMGetLocalVector(veloc_da,&local_V);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_V);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(veloc_da,Veloc_total,INSERT_VALUES,local_V);
	ierr = DMGlobalToLocalEnd(  veloc_da,Veloc_total,INSERT_VALUES,local_V);
	
	ierr = DMDAVecGetArray(veloc_da,local_V,&VV);CHKERRQ(ierr);
	
	
	
	////////
	
	PetscScalar					***TTC;
	Vec						local_TC;
	
	ierr = DMGetLocalVector(thermal_da,&local_TC);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_TC);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(thermal_da,Temper_Cond,INSERT_VALUES,local_TC);
	ierr = DMGlobalToLocalEnd(  thermal_da,Temper_Cond,INSERT_VALUES,local_TC);
	
	ierr = DMDAVecGetArray(thermal_da,local_TC,&TTC);CHKERRQ(ierr);
	
	////////
	
	
	
	PetscInt               M,N,P;

	
	MatStencil ind1[1],ind[8];
	
	PetscScalar u[8*8],val_cond[1];
	
	PetscFunctionBeginUser;
	ierr = DMDAGetInfo(thermal_da,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	
	
	
	
	
	PetscInt i,j,k,c;
	
	PetscInt       sx,sy,sz,mmx,mmy,mmz;
	
	ierr = DMDAGetCorners(thermal_da,&sx,&sy,&sz,&mmx,&mmy,&mmz);CHKERRQ(ierr);
	
	for (k=sz; k<sz+mmz; k++) {
		for (j=sy; j<sy+mmy; j++) {
			for (i=sx; i<sx+mmx; i++) {
			
				ind1[0].i = i;
				ind1[0].j = j;
				ind1[0].k = k;
				
				val_cond[0] = 1.0-TTC[k][j][i];
				ierr = MatSetValuesStencil(A,1,ind1,1,ind1,val_cond,ADD_VALUES);
				
			}
		}
	}
	
	ierr = DMDAGetElementCorners(thermal_da,&sex,&sey,&sez,&mx,&my,&mz);CHKERRQ(ierr);
	for (ek = sez; ek < sez+mz; ek++) {
		for (ej = sey; ej < sey+my; ej++) {
			for (ei = sex; ei < sex+mx; ei++) {
				
				ind[0].i=ei  ; ind[0].j=ej  ; ind[0].k=ek  ;
				ind[1].i=ei+1; ind[1].j=ej  ; ind[1].k=ek  ;
				ind[2].i=ei  ; ind[2].j=ej+1; ind[2].k=ek  ;
				ind[3].i=ei+1; ind[3].j=ej+1; ind[3].k=ek  ;
				ind[4].i=ei  ; ind[4].j=ej  ; ind[4].k=ek+1;
				ind[5].i=ei+1; ind[5].j=ej  ; ind[5].k=ek+1;
				ind[6].i=ei  ; ind[6].j=ej+1; ind[6].k=ek+1;
				ind[7].i=ei+1; ind[7].j=ej+1; ind[7].k=ek+1;
				
				for (i=0;i<8;i++){
					v_vec_aux_ele[i*3+0] = VV[ind[i].k][ind[i].j][ind[i].i].u;
					v_vec_aux_ele[i*3+1] = VV[ind[i].k][ind[i].j][ind[i].i].v;
					v_vec_aux_ele[i*3+2] = VV[ind[i].k][ind[i].j][ind[i].i].w;
				}
				
				montaKeThermal_simplif(TCe_fut, TKe);//modificar
				
				for (c=0;c<T_NE*T_NE;c++) Ttotal[c] = TMe[c] + alpha_thermal*dt_calor_sec*TCe_fut[c];
				
				
				
				for (i=0;i<8;i++){
					if (TTC[ind[i].k][ind[i].j][ind[i].i]==0){
						for (j=0;j<8;j++) u[i*8+j]=0.0;
					}
					else {
						for (j=0;j<8;j++) u[i*8+j]=Ttotal[i*8+j];
					}
				}
				
				ierr = MatSetValuesStencil(A,8,ind,8,ind,u,ADD_VALUES);
			}
		}
	}
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = DMRestoreLocalVector(veloc_da,&local_V);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

PetscErrorCode AssembleF_Thermal(Vec F,DM thermal_da,PetscReal *TKe,PetscReal *TMe,PetscReal *TFe,
								 DM veloc_da, Vec Veloc_total)
{
	//DM                     cda;
	
	Vec                    local_F,local_Temper;
	PetscScalar              ***ff,***tt;
	PetscInt               M,N,P;
	PetscErrorCode         ierr;
	
	PetscInt               sex,sey,sez,mx,my,mz;
	PetscInt               ei,ej,ek;
	
	
	PetscFunctionBeginUser;
	ierr = DMDAGetInfo(thermal_da,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	
	//ierr = DMGetCoordinateDM(thermal_da,&cda);CHKERRQ(ierr);
	
	
	Stokes					***VV;
	Vec						local_V;
	
	ierr = DMGetLocalVector(veloc_da,&local_V);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_V);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(veloc_da,Veloc_total,INSERT_VALUES,local_V);
	ierr = DMGlobalToLocalEnd(  veloc_da,Veloc_total,INSERT_VALUES,local_V);
	
	ierr = DMDAVecGetArray(veloc_da,local_V,&VV);CHKERRQ(ierr);
	
	
	////////
	
	PetscScalar					***TTC;
	Vec						local_TC;
	
	ierr = DMGetLocalVector(thermal_da,&local_TC);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_TC);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(thermal_da,Temper_Cond,INSERT_VALUES,local_TC);
	ierr = DMGlobalToLocalEnd(  thermal_da,Temper_Cond,INSERT_VALUES,local_TC);
	
	ierr = DMDAVecGetArray(thermal_da,local_TC,&TTC);CHKERRQ(ierr);
	
	////////
	
	
	
	/* get acces to the vector */
	ierr = DMGetLocalVector(thermal_da,&local_F);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_F);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(thermal_da,local_F,&ff);CHKERRQ(ierr);
	
	
	ierr = DMGetLocalVector(thermal_da,&local_Temper);CHKERRQ(ierr);
	
	ierr = VecZeroEntries(local_Temper);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(thermal_da,Temper,INSERT_VALUES,local_Temper);
	ierr = DMGlobalToLocalEnd(thermal_da,Temper,INSERT_VALUES,local_Temper);
	
	ierr = DMDAVecGetArray(thermal_da,local_Temper,&tt);CHKERRQ(ierr);
	
	PetscInt i,j,k,c;
	
	MatStencil ind[8];
	
	
	ierr = DMDAGetElementCorners(thermal_da,&sex,&sey,&sez,&mx,&my,&mz);CHKERRQ(ierr);
	for (ek = sez; ek < sez+mz; ek++) {
		for (ej = sey; ej < sey+my; ej++) {
			for (ei = sex; ei < sex+mx; ei++) {
				
				
				ind[0].i=ei  ; ind[0].j=ej  ; ind[0].k=ek  ;
				ind[1].i=ei+1; ind[1].j=ej  ; ind[1].k=ek  ;
				ind[2].i=ei  ; ind[2].j=ej+1; ind[2].k=ek  ;
				ind[3].i=ei+1; ind[3].j=ej+1; ind[3].k=ek  ;
				ind[4].i=ei  ; ind[4].j=ej  ; ind[4].k=ek+1;
				ind[5].i=ei+1; ind[5].j=ej  ; ind[5].k=ek+1;
				ind[6].i=ei  ; ind[6].j=ej+1; ind[6].k=ek+1;
				ind[7].i=ei+1; ind[7].j=ej+1; ind[7].k=ek+1;
				
				
				
				for (i=0;i<8;i++){
					v_vec_aux_ele[i*3+0] = VV[ind[i].k][ind[i].j][ind[i].i].u;
					v_vec_aux_ele[i*3+1] = VV[ind[i].k][ind[i].j][ind[i].i].v;
					v_vec_aux_ele[i*3+2] = VV[ind[i].k][ind[i].j][ind[i].i].w;
				}
				
				
				montaKeThermal_simplif(TCe, TKe);//modificar
				
				for (c=0;c<T_NE*T_NE;c++) Ttotal_b[c] = TMe[c] - comp_alpha_thermal*dt_calor_sec*TCe[c];
				
				/*if (ek==sez && ej==sey && ei==sex){
					printf("Ttotal_b\n");
					for (c=0;c<T_NE;c++){
						for (j=0;j<T_NE;j++){
							printf("%5.1g ",Ttotal_b[c*T_NE+j]);
						}
						printf("\n");
					}
					printf("\n");
				}*/
				
				
				
				
				///VecGetValues(T_vec,T_NE,Indices_Ke_Thermal_I,T_vec_aux_ele);
				
				//for (c=0;c<T_NE;c++)
				//	printf("%5.1lg ",tt[ind[c].i][ind[c].j][ind[c].k]);
				//printf("\n");
				
				for (c=0;c<T_NE;c++){
					T_vec_aux_ele_final[c]=dt_calor_sec*TFe[c];
					for (j=0;j<T_NE;j++){
						T_vec_aux_ele_final[c]+=tt[ind[j].k][ind[j].j][ind[j].i]*Ttotal_b[c*T_NE+j];
					}
				}
				
				for (i=0;i<8;i++){
					//if (ind[i].k==0 || ind[i].k==P-1){
					if (TTC[ind[i].k][ind[i].j][ind[i].i]==0){
						T_vec_aux_ele_final[i]=0.0;
					}
				}
				
				for (c=0;c<T_NE;c++){
					ff[ind[c].k][ind[c].j][ind[c].i] += T_vec_aux_ele_final[c];
				}
				
				
				
				/*for (i=0;i<8;i++){
					if (ind[i].k==0 || ind[i].k==P-1){
						for (j=0;j<8;j++) u[i*8+j]=0.0;
					}
					else {
						for (j=0;j<8;j++) u[i*8+j]=TKe[i*8+j];
					}
				}*/
				
				
			}
		}
	}
	
	
	
	PetscInt       sx,sy,sz,mmx,mmy,mmz;
	
	ierr = DMDAGetCorners(thermal_da,&sx,&sy,&sz,&mmx,&mmy,&mmz);CHKERRQ(ierr);
	
	for (k=sz; k<sz+mmz; k++) {
		for (j=sy; j<sy+mmy; j++) {
			for (i=sx; i<sx+mmx; i++) {
				//if (k==0 || k==P-1){
				if (TTC[k][j][i]==0){
					ff[k][j][i]=tt[k][j][i];
				}
			}
		}
	}
	
	ierr = DMDAVecRestoreArray(thermal_da,local_F,&ff);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(thermal_da,local_F,ADD_VALUES,F);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(thermal_da,local_F,ADD_VALUES,F);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(thermal_da,&local_F);CHKERRQ(ierr);
	
	ierr = DMRestoreLocalVector(thermal_da,&local_Temper);CHKERRQ(ierr);
	
	
	ierr = DMRestoreLocalVector(veloc_da,&local_V);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#include <math.h>


PetscErrorCode Thermal_init(Vec F,DM thermal_da)
{
	//DM                     cda;
	
	Vec                    local_F;
	PetscScalar              ***ff;
	PetscInt               M,N,P;
	PetscErrorCode         ierr;
	
	PetscInt i,j,k;
	
	
	PetscFunctionBeginUser;
	ierr = DMDAGetInfo(thermal_da,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	
	//ierr = DMGetCoordinateDM(thermal_da,&cda);CHKERRQ(ierr);
	
	
	
	/* get acces to the vector */
	ierr = DMGetLocalVector(thermal_da,&local_F);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_F);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(thermal_da,local_F,&ff);CHKERRQ(ierr);
	
	PetscInt       sx,sy,sz,mmx,mmy,mmz;
	
	ierr = DMDAGetCorners(thermal_da,&sx,&sy,&sz,&mmx,&mmy,&mmz);CHKERRQ(ierr);
	
	PetscReal temper_aux,t1_aux;
	
	for (k=sz; k<sz+mmz; k++) {
		for (j=sy; j<sy+mmy; j++) {
			for (i=sx; i<sx+mmx; i++) {
				
				
				if (T_initial_cond==0){
					temper_aux=(Delta_T*(P-1-k))/(P-1) + 100*cos(i*3.14159/(M-1));
				}
				
				if (T_initial_cond==1){
					temper_aux=(Delta_T*(P-1-k))/(P-1) + 100*cos(j*3.14159/(N-1))*cos(i*3.14159/(M-1));
				}
				
				if (T_initial_cond==74){
					t1_aux = 0.5*tan((P-1-k)*3./(P-1)-1.5)/tan(1.5)+0.5;
					temper_aux=Delta_T*t1_aux + Delta_T*0.2*tan(i*3./(M-1)-1.5);
				}
				
				
				if (k==P-1) temper_aux=0.0;
				if (k==0) temper_aux=Delta_T;
				
				if (temper_aux>Delta_T) temper_aux=Delta_T;
				if (temper_aux<0.0) temper_aux=0.0;
				
				ff[k][j][i]=temper_aux;
			}
		}
	}
	
	
	
	ierr = DMDAVecRestoreArray(thermal_da,local_F,&ff);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(thermal_da,local_F,ADD_VALUES,F);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(thermal_da,local_F,ADD_VALUES,F);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(thermal_da,&local_F);CHKERRQ(ierr);
	
	//ierr = DMDAVecRestoreArray(cda,coords,&_coords);CHKERRQ(ierr);
	PetscFunctionReturn(0);
	
}

