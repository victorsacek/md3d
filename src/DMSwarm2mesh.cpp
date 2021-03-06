#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmswarm.h>

extern DM dms;

extern DM da_Thermal;

extern Vec geoq,local_geoq;
extern Vec geoq_rho,local_geoq_rho;
extern Vec geoq_cont,local_geoq_cont;
extern Vec geoq_H,local_geoq_H;
extern Vec geoq_strain, local_geoq_strain;
extern Vec geoq_strain_rate, local_geoq_strain_rate;

extern Vec Temper,local_Temper;

extern long Nx,Ny,Nz;

extern double dx_const;
extern double dy_const;
extern double dz_const;

extern double Lx, Ly, depth;

extern PetscInt visc_harmonic_mean;
extern PetscInt visc_const_per_element;


PetscErrorCode Swarm2Mesh(){

	PetscErrorCode ierr;
	PetscScalar             ***qq,***qq_cont,***qq_rho,***TT,***qq_H,***qq_strain;
	PetscScalar				***qq_strain_rate;
	
	ierr = VecSet(geoq,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_rho,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_cont,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_H,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_strain,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_strain_rate,0.0);CHKERRQ(ierr);
	
	ierr = VecZeroEntries(local_geoq);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_geoq_rho);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_geoq_H);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_geoq_strain);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_geoq_strain_rate);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(da_Thermal,geoq,INSERT_VALUES,local_geoq);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq,INSERT_VALUES,local_geoq);
	
	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
	
	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_H,INSERT_VALUES,local_geoq_H);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_H,INSERT_VALUES,local_geoq_H);
	
	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_strain,INSERT_VALUES,local_geoq_strain);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_strain,INSERT_VALUES,local_geoq_strain);

	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_strain_rate,INSERT_VALUES,local_geoq_strain_rate);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_strain_rate,INSERT_VALUES,local_geoq_strain_rate);
	
	
	ierr = DMDAVecGetArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_H,&qq_H);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_strain,&qq_strain);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_strain_rate,&qq_strain_rate);CHKERRQ(ierr);

	
	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);
	
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
	
	

	PetscInt nlocal,bs,p;
	
	PetscReal *array;
	PetscReal *geoq_fac;
	PetscReal *rho_fac;
	PetscReal *H_fac;
	PetscReal *strain_fac;
	PetscReal *strain_rate_fac;
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	
	ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"geoq_fac",NULL,NULL,(void**)&geoq_fac);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"rho_fac",NULL,NULL,(void**)&rho_fac);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"H_fac",NULL,NULL,(void**)&H_fac);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"strain_fac",NULL,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"strain_rate_fac",NULL,NULL,(void**)&strain_rate_fac);CHKERRQ(ierr);
	
	if (visc_harmonic_mean==1){
		for (p=0; p<nlocal; p++) {
			PetscReal cx,cy,cz;
			PetscReal rx,ry,rz,rfac;
			PetscInt i,j,k;
			
			cx = array[3*p];
			cy = array[3*p+1];
			cz = array[3*p+2];
			
			i = (int)(cx/dx_const);
			j = (int)(cy/dy_const);
			k = (int)((cz+depth)/dz_const);
			
			
			
			if (i<0 || i>=Nx-1) {printf("estranho i=%d\n",i); exit(1);}
			if (j<0 || j>=Ny-1) {printf("estranho j=%d\n",j); exit(1);}
			if (k<0 || k>=Nz-1) {printf("estranho k=%d\n",k); exit(1);}
			
			if (i==Nx-1) i=Nx-2;
			if (j==Ny-1) j=Ny-2;
			if (k==Nz-1) k=Nz-2;
			
			rx = (cx-i*dx_const)/dx_const;
			ry = (cy-j*dy_const)/dy_const;
			rz = (cz-(-depth+k*dz_const))/dz_const;
			
			if (rx<0 || rx>1) {printf("estranho rx=%f\n",rx); exit(1);}
			if (ry<0 || ry>1) {printf("estranho ry=%f\n",ry); exit(1);}
			if (rz<0 || rz>1) {printf("estranho rz=%f\n",rz); exit(1);}
			
			
			rfac = (1.0-rx)*(1.0-ry)*(1.0-rz);
			qq		[k][j][i] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k][j][i] += rfac*rho_fac[p];
			qq_H	[k][j][i] += rfac*H_fac[p];
			qq_strain[k][j][i] += rfac*strain_fac[p];
			qq_strain_rate[k][j][i] += rfac*strain_rate_fac[p];
			qq_cont	[k][j][i] += rfac;
			
			rfac = (rx)*(1.0-ry)*(1.0-rz);
			qq		[k][j][i+1] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k][j][i+1] += rfac*rho_fac[p];
			qq_H	[k][j][i+1] += rfac*H_fac[p];
			qq_strain[k][j][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k][j][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k][j][i+1] += rfac;
			
			rfac = (1.0-rx)*(ry)*(1.0-rz);
			qq		[k][j+1][i] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k][j+1][i] += rfac*rho_fac[p];
			qq_H	[k][j+1][i] += rfac*H_fac[p];
			qq_strain[k][j+1][i] += rfac*strain_fac[p];
			qq_strain_rate[k][j+1][i] += rfac*strain_rate_fac[p];
			qq_cont	[k][j+1][i] += rfac;
			
			rfac = (rx)*(ry)*(1.0-rz);
			qq		[k][j+1][i+1] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k][j+1][i+1] += rfac*rho_fac[p];
			qq_H	[k][j+1][i+1] += rfac*H_fac[p];
			qq_strain[k][j+1][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k][j+1][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k][j+1][i+1] += rfac;
			
			
			
			rfac = (1.0-rx)*(1.0-ry)*(rz);
			qq		[k+1][j][i] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k+1][j][i] += rfac*rho_fac[p];
			qq_H	[k+1][j][i] += rfac*H_fac[p];
			qq_strain[k+1][j][i] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j][i] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j][i] += rfac;
			
			rfac = (rx)*(1.0-ry)*(rz);
			qq		[k+1][j][i+1] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k+1][j][i+1] += rfac*rho_fac[p];
			qq_H	[k+1][j][i+1] += rfac*H_fac[p];
			qq_strain[k+1][j][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j][i+1] += rfac;
			
			rfac = (1.0-rx)*(ry)*(rz);
			qq		[k+1][j+1][i] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k+1][j+1][i] += rfac*rho_fac[p];
			qq_H	[k+1][j+1][i] += rfac*H_fac[p];
			qq_strain[k+1][j+1][i] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j+1][i] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j+1][i] += rfac;
			
			rfac = (rx)*(ry)*(rz);
			qq		[k+1][j+1][i+1] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k+1][j+1][i+1] += rfac*rho_fac[p];
			qq_H	[k+1][j+1][i+1] += rfac*H_fac[p];
			qq_strain[k+1][j+1][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j+1][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j+1][i+1] += rfac;
		}
	}
	else{
		for (p=0; p<nlocal; p++) {
			PetscReal cx,cy,cz;
			PetscReal rx,ry,rz,rfac;
			PetscInt i,j,k;
			
			cx = array[3*p];
			cy = array[3*p+1];
			cz = array[3*p+2];
			
			i = (int)(cx/dx_const);
			j = (int)(cy/dy_const);
			k = (int)((cz+depth)/dz_const);
			
			
			
			if (i<0 || i>=Nx-1) {printf("estranho i=%d\n",i); exit(1);}
			if (j<0 || j>=Ny-1) {printf("estranho j=%d\n",j); exit(1);}
			if (k<0 || k>=Nz-1) {printf("estranho k=%d\n",k); exit(1);}
			
			if (i==Nx-1) i=Nx-2;
			if (j==Ny-1) j=Ny-2;
			if (k==Nz-1) k=Nz-2;
			
			rx = (cx-i*dx_const)/dx_const;
			ry = (cy-j*dy_const)/dy_const;
			rz = (cz-(-depth+k*dz_const))/dz_const;
			
			if (rx<0 || rx>1) {printf("estranho rx=%f\n",rx); exit(1);}
			if (ry<0 || ry>1) {printf("estranho ry=%f\n",ry); exit(1);}
			if (rz<0 || rz>1) {printf("estranho rz=%f\n",rz); exit(1);}
			
			
			rfac = (1.0-rx)*(1.0-ry)*(1.0-rz);
			qq		[k][j][i] += rfac*geoq_fac[p];
			qq_rho	[k][j][i] += rfac*rho_fac[p];
			qq_H	[k][j][i] += rfac*H_fac[p];
			qq_strain[k][j][i] += rfac*strain_fac[p];
			qq_strain_rate[k][j][i] += rfac*strain_rate_fac[p];
			qq_cont	[k][j][i] += rfac;
			
			rfac = (rx)*(1.0-ry)*(1.0-rz);
			qq		[k][j][i+1] += rfac*geoq_fac[p];
			qq_rho	[k][j][i+1] += rfac*rho_fac[p];
			qq_H	[k][j][i+1] += rfac*H_fac[p];
			qq_strain[k][j][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k][j][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k][j][i+1] += rfac;
			
			rfac = (1.0-rx)*(ry)*(1.0-rz);
			qq		[k][j+1][i] += rfac*geoq_fac[p];
			qq_rho	[k][j+1][i] += rfac*rho_fac[p];
			qq_H	[k][j+1][i] += rfac*H_fac[p];
			qq_strain[k][j+1][i] += rfac*strain_fac[p];
			qq_strain_rate[k][j+1][i] += rfac*strain_rate_fac[p];
			qq_cont	[k][j+1][i] += rfac;
			
			rfac = (rx)*(ry)*(1.0-rz);
			qq		[k][j+1][i+1] += rfac*geoq_fac[p];
			qq_rho	[k][j+1][i+1] += rfac*rho_fac[p];
			qq_H	[k][j+1][i+1] += rfac*H_fac[p];
			qq_strain[k][j+1][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k][j+1][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k][j+1][i+1] += rfac;
			
			
			
			rfac = (1.0-rx)*(1.0-ry)*(rz);
			qq		[k+1][j][i] += rfac*geoq_fac[p];
			qq_rho	[k+1][j][i] += rfac*rho_fac[p];
			qq_H	[k+1][j][i] += rfac*H_fac[p];
			qq_strain[k+1][j][i] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j][i] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j][i] += rfac;
			
			rfac = (rx)*(1.0-ry)*(rz);
			qq		[k+1][j][i+1] += rfac*geoq_fac[p];
			qq_rho	[k+1][j][i+1] += rfac*rho_fac[p];
			qq_H	[k+1][j][i+1] += rfac*H_fac[p];
			qq_strain[k+1][j][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j][i+1] += rfac;
			
			rfac = (1.0-rx)*(ry)*(rz);
			qq		[k+1][j+1][i] += rfac*geoq_fac[p];
			qq_rho	[k+1][j+1][i] += rfac*rho_fac[p];
			qq_H	[k+1][j+1][i] += rfac*H_fac[p];
			qq_strain[k+1][j+1][i] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j+1][i] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j+1][i] += rfac;
			
			rfac = (rx)*(ry)*(rz);
			qq		[k+1][j+1][i+1] += rfac*geoq_fac[p];
			qq_rho	[k+1][j+1][i+1] += rfac*rho_fac[p];
			qq_H	[k+1][j+1][i+1] += rfac*H_fac[p];
			qq_strain[k+1][j+1][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j+1][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j+1][i+1] += rfac;
		}
	}
	
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq,ADD_VALUES,geoq);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq,ADD_VALUES,geoq);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_rho,ADD_VALUES,geoq_rho);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_rho,ADD_VALUES,geoq_rho);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_H,&qq_H);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_H,ADD_VALUES,geoq_H);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_H,ADD_VALUES,geoq_H);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_strain,&qq_strain);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_strain,ADD_VALUES,geoq_strain);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_strain,ADD_VALUES,geoq_strain);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_strain_rate,&qq_strain_rate);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_strain_rate,ADD_VALUES,geoq_strain_rate);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_strain_rate,ADD_VALUES,geoq_strain_rate);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_cont,ADD_VALUES,geoq_cont);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_cont,ADD_VALUES,geoq_cont);CHKERRQ(ierr);
	
	//VecPointwiseMax(geoq_cont,geoq_cont,geoqOnes);
	VecPointwiseDivide(geoq,geoq,geoq_cont);
	if (visc_harmonic_mean==1) VecReciprocal(geoq); //<--- harmonic
		
	VecPointwiseDivide(geoq_rho,geoq_rho,geoq_cont);
	VecPointwiseDivide(geoq_H,geoq_H,geoq_cont);
	VecPointwiseDivide(geoq_strain,geoq_strain,geoq_cont);
	VecPointwiseDivide(geoq_strain_rate,geoq_strain_rate,geoq_cont);
	//VecPointwiseMax(geoq,geoq,geoqOnes);


	if (visc_const_per_element==1){
		ierr = VecSet(geoq,0.0);CHKERRQ(ierr);
		ierr = VecSet(geoq_cont,0.0);CHKERRQ(ierr);

		ierr = DMGlobalToLocalBegin(da_Thermal,geoq,INSERT_VALUES,local_geoq);
		ierr = DMGlobalToLocalEnd(  da_Thermal,geoq,INSERT_VALUES,local_geoq);
		ierr = DMDAVecGetArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);

		
		ierr = DMGlobalToLocalBegin(da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);
		ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);
		ierr = DMDAVecGetArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
		if (visc_harmonic_mean==1){
			for (p=0; p<nlocal; p++) {
				PetscReal cx,cy,cz;
				PetscInt i,j,k;
				
				cx = array[3*p];
				cy = array[3*p+1];
				cz = array[3*p+2];
				
				i = (int)(cx/dx_const);
				j = (int)(cy/dy_const);
				k = (int)((cz+depth)/dz_const);
				
				qq[k][j][i] += 1.0/geoq_fac[p];
				qq_cont[k][j][i] += 1.0;
			}
		}
		else {
			for (p=0; p<nlocal; p++) {
				PetscReal cx,cy,cz;
				PetscInt i,j,k;
				
				cx = array[3*p];
				cy = array[3*p+1];
				cz = array[3*p+2];
				
				i = (int)(cx/dx_const);
				j = (int)(cy/dy_const);
				k = (int)((cz+depth)/dz_const);
				
				qq[k][j][i] += geoq_fac[p];
				qq_cont[k][j][i] += 1.0;
			}
		}


		ierr = DMDAVecRestoreArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq,ADD_VALUES,geoq);CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq,ADD_VALUES,geoq);CHKERRQ(ierr);
		
		ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_cont,ADD_VALUES,geoq_cont);CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_cont,ADD_VALUES,geoq_cont);CHKERRQ(ierr);

		VecPointwiseDivide(geoq,geoq,geoq_cont);
		if (visc_harmonic_mean==1) VecReciprocal(geoq); //<--- harmonic

	}

	
	


	ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"geoq_fac",NULL,NULL,(void**)&geoq_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"rho_fac",NULL,NULL,(void**)&rho_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"H_fac",NULL,NULL,(void**)&H_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_fac",NULL,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_rate_fac",NULL,NULL,(void**)&strain_rate_fac);CHKERRQ(ierr);

	ierr = DMGlobalToLocalBegin(da_Thermal,Temper,INSERT_VALUES,local_Temper);
	ierr = DMGlobalToLocalEnd(  da_Thermal,Temper,INSERT_VALUES,local_Temper);
	
	ierr = DMDAVecGetArray(da_Thermal,local_Temper,&TT);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
	
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);

	
	PetscInt       sx,sy,sz,mmx,mmy,mmz;
	
	ierr = DMDAGetCorners(da_Thermal,&sx,&sy,&sz,&mmx,&mmy,&mmz);CHKERRQ(ierr);
	
	int k,j,i;
	
	for (k=sz; k<sz+mmz; k++) {
		for (j=sy; j<sy+mmy; j++) {
			for (i=sx; i<sx+mmx; i++) {
				if (qq_rho[k][j][i]<100.0){
					TT[k][j][i]=0.0;
				}
			}
		}
	}
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_Temper,&TT);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_Temper,INSERT_VALUES,Temper);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_Temper,INSERT_VALUES,Temper);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
	
}
