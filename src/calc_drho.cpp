#include <petscksp.h>
#include <petscdmda.h>

extern Vec dRho;
extern Vec Temper;
extern Vec geoq_rho;

extern double alpha_exp_thermo;
extern double gravity;
extern double RHOM;

PetscErrorCode calc_drho()
{
	PetscErrorCode ierr=0;
	
	VecCopy(Temper, dRho); CHKERRQ(ierr);
	
	PetscScalar a;
	
	a=alpha_exp_thermo*gravity*RHOM; CHKERRQ(ierr);
	
	//VecScale(dRho,a); CHKERRQ(ierr);
	VecAXPBY(dRho,-gravity,a,geoq_rho); CHKERRQ(ierr);
	
	return ierr;
	
}

PetscErrorCode write_pressure(int cont);

extern Vec local_P;
extern Vec local_P_aux;
extern Vec Pressure;

extern Vec Pressure_aux;

extern Vec local_geoq_rho;

extern DM da_Veloc;

extern DM da_Thermal;

extern long Nx,Ny,Nz;

extern long V_GT;
extern long T_NE;

extern double Lx, Ly, depth;


extern int n_interfaces;
extern PetscScalar *interfaces;

extern PetscScalar *inter_rho;

typedef struct {
	PetscScalar u;
	PetscScalar v;
	PetscScalar w;
} Stokes;

PetscErrorCode DMDAGetElementCorners(DM da,PetscInt *sx,PetscInt *sy, PetscInt *sz,PetscInt *mx, PetscInt *my, PetscInt *mz);

PetscErrorCode calc_pressure()
{
	Stokes					***pp;
	
	PetscErrorCode         ierr;
	
	ierr = VecZeroEntries(local_P);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(da_Veloc,Pressure,INSERT_VALUES,local_P);
	ierr = DMGlobalToLocalEnd(  da_Veloc,Pressure,INSERT_VALUES,local_P);
	
	ierr = DMDAVecGetArray(da_Veloc,local_P,&pp);CHKERRQ(ierr);
	
	
	
	//PetscScalar             **qq_rho;
	
	//ierr = DMGlobalToLocalBegin(da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
	//ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
	
	//ierr = DMDAVecGetArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);
	
	PetscInt i,j,k,c,n,g;
	
	MatStencil indr[T_NE],ind[V_GT];
	
	int M,N,P;
	double xx,yy,zz;
	
	PetscFunctionBeginUser;
	ierr = DMDAGetInfo(da_Thermal,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	
	printf("%d %d %d P\n",M,N,P);
	
	//
	
	PetscInt               sex1,sey1,sez1,mx1,my1,mz1;
	
	PetscInt               sex,sey,sez,mx,my,mz;
	PetscInt               ei,ej,ek;
	
	ierr = DMDAGetElementCorners(da_Thermal,&sex,&sey,&sez,&mx,&my,&mz);CHKERRQ(ierr);
	
	ierr = DMDAGetElementCorners(da_Veloc,&sex1,&sey1,&sez1,&mx1,&my1,&mz1);CHKERRQ(ierr);
	
	if ((sex-sex1!=0)||(sey-sey1!=0)||(sez-sez1!=0)||(mx-mx1!=0)||(my-my1!=0)||(mz-mz1!=0)){
		printf("%d %d %d %d %d %d\n",sex-sex1,sey-sey1,sez-sez1,mx-mx1,my-my1,mz-mz1);
		SETERRQ1(PETSC_COMM_WORLD,1,"Wrong partition (temper,velocity)\n",1);
	}
	
	PetscReal interp_interfaces[n_interfaces];
	
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
				
				xx = ei*Lx/(M-1);
				yy = ej*Ly/(N-1);
				zz = -(P-1-ek)*depth/(P-1);
				
				for (int in=0;in<n_interfaces;in++){
					float rfac = 0.25;
					interp_interfaces[in] = interfaces[ej*Nx+ei + Nx*Ny*in] * rfac;
					
					rfac = 0.25;
					interp_interfaces[in] += interfaces[ej*Nx+(ei+1) + Nx*Ny*in] * rfac;
					
					rfac = 0.25;
					interp_interfaces[in] += interfaces[(ej+1)*Nx+ei + Nx*Ny*in] * rfac;
					
					rfac = 0.25;
					interp_interfaces[in] += interfaces[(ej+1)*Nx+(ei+1) + Nx*Ny*in] * rfac;
				}
				
				double pressure_cumulat = 0;
				
				int verif=0;
				
				
				
				if (zz>interp_interfaces[n_interfaces-1]){
					pressure_cumulat += inter_rho[n_interfaces]*gravity*(-zz);
				}
				else{
					pressure_cumulat += inter_rho[n_interfaces]*gravity*
										(-interp_interfaces[n_interfaces-1]);
					for (int in=n_interfaces-1; in>0 && verif==0;in--){
						if (zz>interp_interfaces[in-1]){
							pressure_cumulat += inter_rho[in]*gravity*
												(interp_interfaces[in]-zz);
							verif=1;
						}
						else{
							pressure_cumulat += inter_rho[in]*gravity*
												(interp_interfaces[in]-interp_interfaces[in-1]);
						}
					}
					if (verif==0){
						pressure_cumulat += inter_rho[0]*gravity*(interp_interfaces[0]-zz);
					}
				}
				
				pp[ek][ej][ei].u = pressure_cumulat;
				
				//printf("%d ",in);
				//if (ek==1) printf("pp = %lf\n",pressure_cumulat);
				
			}
		}
	}
	
	
	//ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Veloc,local_P,&pp);CHKERRQ(ierr);
	
	ierr = DMLocalToGlobalBegin(da_Veloc,local_P,INSERT_VALUES,Pressure);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Veloc,local_P,INSERT_VALUES,Pressure);CHKERRQ(ierr);
	
	//write_pressure(-1);
	
	return ierr;
	
}



PetscErrorCode shift_pressure() //necessary if the surface pressure is not close to zero
{
	PetscErrorCode         ierr;

	Stokes					***pp;
	PetscScalar				***pp_aux;

	//get Pressure array
	ierr = DMGlobalToLocalBegin(da_Veloc,Pressure,INSERT_VALUES,local_P);
	ierr = DMGlobalToLocalEnd(  da_Veloc,Pressure,INSERT_VALUES,local_P);
	ierr = DMDAVecGetArray(da_Veloc,local_P,&pp);CHKERRQ(ierr);

	//get Pressure_aux array
	ierr = DMGlobalToLocalBegin(da_Thermal,Pressure_aux,INSERT_VALUES,local_P_aux);
	ierr = DMGlobalToLocalEnd(  da_Thermal,Pressure_aux,INSERT_VALUES,local_P_aux);
	ierr = DMDAVecGetArray(da_Thermal,local_P_aux,&pp_aux);CHKERRQ(ierr);

	
	PetscInt       sx,sy,sz,mmx,mmy,mmz;
	PetscInt i,j,k;
	
	ierr = DMDAGetCorners(da_Veloc,&sx,&sy,&sz,&mmx,&mmy,&mmz);CHKERRQ(ierr);

	PetscScalar ppp;
	PetscInt cont_ppp;

	for (k=sz; k<sz+mmz; k++) {
		for (j=sy; j<sy+mmy; j++) {
			for (i=sx; i<sx+mmx; i++) {
				ppp=0.0;
				cont_ppp = 0;
				if (k<Nz-1 && j<Ny-1 && i<Nx-1) {ppp+=pp[k  ][j  ][i  ].u; cont_ppp++;}
				if (k<Nz-1 && j<Ny-1 && i>0   )	{ppp+=pp[k  ][j  ][i-1].u; cont_ppp++;}
				if (k<Nz-1 && j>0    && i<Nx-1) {ppp+=pp[k  ][j-1][i  ].u; cont_ppp++;}
				if (k<Nz-1 && j>0    && i>0   )	{ppp+=pp[k  ][j-1][i-1].u; cont_ppp++;}
				if (k>0    && j<Ny-1 && i<Nx-1)	{ppp+=pp[k-1][j  ][i  ].u; cont_ppp++;}
				if (k>0    && j<Ny-1 && i>0   )	{ppp+=pp[k-1][j  ][i-1].u; cont_ppp++;}
				if (k>0    && j>0    && i<Nx-1)	{ppp+=pp[k-1][j-1][i  ].u; cont_ppp++;}
				if (k>0    && j>0    && i>0   )	{ppp+=pp[k-1][j-1][i-1].u; cont_ppp++;}
				pp_aux[k][j][i] = ppp/cont_ppp;
			}
		}
	}
	int cont_air = 0, cont_air_all;
	float pressure_air = 0.0, pressure_air_all;

	for (k=sz; k<sz+mmz; k++) {
		for (j=sy; j<sy+mmy; j++) {
			for (i=sx; i<sx+mmx; i++) {
				if (k==Nz-1){
					pressure_air+=pp_aux[k][j][i];
					cont_air++;
				}
			}
		}
	}

	//restore Pressure_aux
	ierr = DMDAVecRestoreArray(da_Thermal,local_P_aux,&pp_aux);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_P_aux,INSERT_VALUES,Pressure_aux);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_P_aux,INSERT_VALUES,Pressure_aux);CHKERRQ(ierr);

	MPI_Allreduce(&cont_air,&cont_air_all,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
	MPI_Allreduce(&pressure_air,&pressure_air_all,1,MPI_FLOAT,MPI_SUM,PETSC_COMM_WORLD);

	PetscReal pressure_min;
	
	pressure_min = pressure_air_all/cont_air_all;
	pressure_min=-pressure_min;
	VecShift(Pressure_aux,pressure_min);


	//restore Pressure
	ierr = DMDAVecRestoreArray(da_Veloc,local_P,&pp);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Veloc,local_P,INSERT_VALUES,Pressure);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Veloc,local_P,INSERT_VALUES,Pressure);CHKERRQ(ierr);



	return ierr;
	
}