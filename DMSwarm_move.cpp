#include <petscdmda.h>
#include <petscdmswarm.h>

typedef struct {
	PetscScalar u;
	PetscScalar v;
	PetscScalar w;
	//PetscScalar p;
} Stokes;

//PetscErrorCode SwarmViewGP(DM dms,const char prefix[]);

extern DM dms;

extern DM da_Veloc;

extern long Nx,Ny,Nz;

extern double dx_const;
extern double dy_const;
extern double dz_const;

extern double Lx, Ly, depth;

extern Vec local_V,Veloc_fut;

PetscErrorCode moveSwarm(PetscReal dt)
{
	PetscErrorCode ierr=0;
	
	Stokes					***VV;
	ierr = VecZeroEntries(local_V);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(da_Veloc,Veloc_fut,INSERT_VALUES,local_V);
	ierr = DMGlobalToLocalEnd(  da_Veloc,Veloc_fut,INSERT_VALUES,local_V);
	
	ierr = DMDAVecGetArray(da_Veloc,local_V,&VV);CHKERRQ(ierr);
	
	
	
	PetscInt nlocal,bs,p;
	
	PetscReal *array;
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	
	ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	
	for (p=0; p<nlocal; p++) {
		PetscReal cx,cy,cz,vx,vy,vz;
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
		
		
		
		//VV[k][j][i].u + ;
		
		rx = (cx-i*dx_const)/dx_const;
		ry = (cy-j*dy_const)/dy_const;
		rz = (cz-(-depth+k*dz_const))/dz_const;
		
		if (rx<0 || rx>1) {printf("estranho rx=%f\n",rx); exit(1);}
		if (ry<0 || ry>1) {printf("estranho ry=%f\n",ry); exit(1);}
		if (rz<0 || rz>1) {printf("estranho rz=%f\n",rz); exit(1);}
		
		rfac = (1.0-rx)*(1.0-ry)*(1.0-rz);
		vx = VV[k][j][i].u * rfac;
		vy = VV[k][j][i].v * rfac;
		vz = VV[k][j][i].w * rfac;
		
		rfac = (rx)*(1.0-ry)*(1.0-rz);
		vx += VV[k][j][i+1].u * rfac;
		vy += VV[k][j][i+1].v * rfac;
		vz += VV[k][j][i+1].w * rfac;
		
		rfac = (1.0-rx)*(ry)*(1.0-rz);
		vx += VV[k][j+1][i].u * rfac;
		vy += VV[k][j+1][i].v * rfac;
		vz += VV[k][j+1][i].w * rfac;
		
		rfac = (rx)*(ry)*(1.0-rz);
		vx += VV[k][j+1][i+1].u * rfac;
		vy += VV[k][j+1][i+1].v * rfac;
		vz += VV[k][j+1][i+1].w * rfac;
		
		rfac = (1.0-rx)*(1.0-ry)*(rz);
		vx += VV[k+1][j][i].u * rfac;
		vy += VV[k+1][j][i].v * rfac;
		vz += VV[k+1][j][i].w * rfac;
		
		rfac = (rx)*(1.0-ry)*(rz);
		vx += VV[k+1][j][i+1].u * rfac;
		vy += VV[k+1][j][i+1].v * rfac;
		vz += VV[k+1][j][i+1].w * rfac;
		
		rfac = (1.0-rx)*(ry)*(rz);
		vx += VV[k+1][j+1][i].u * rfac;
		vy += VV[k+1][j+1][i].v * rfac;
		vz += VV[k+1][j+1][i].w * rfac;
		
		rfac = (rx)*(ry)*(rz);
		vx += VV[k+1][j+1][i+1].u * rfac;
		vy += VV[k+1][j+1][i+1].v * rfac;
		vz += VV[k+1][j+1][i+1].w * rfac;
		
		
		
		//dt=0.01;
		//vx = 0;
		//vy = -cz;
		//vz = cy;
		
		array[3*p  ] += dt * vx;
		array[3*p+1] += dt * vy;
		array[3*p+2] += dt * vz;
		
	}
	ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	
	ierr = DMSwarmMigrate(dms,PETSC_TRUE);CHKERRQ(ierr);
	
	//PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN-1,"step%d",tk);
	//ierr = SwarmViewGP(dms,prefix);CHKERRQ(ierr);
	
	//ierr = SwarmViewGP(dms,"step1");CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Veloc,local_V,&VV);CHKERRQ(ierr);
	
	//exit(1);
	
	PetscFunctionReturn(0);
	
}