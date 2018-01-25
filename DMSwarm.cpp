#include <petscsf.h>
#include <petscdm.h>
#include <petscksp.h>
#include <petscdmda.h>
#include <petscdmshell.h>
#include <petscdmswarm.h>
#include <petsc/private/dmimpl.h>
#include <petscmath.h>
//#include <petscsys.h>

extern DM dmcell;

extern DM dms;

extern DM da_Veloc;

extern double dx_const;
extern double dy_const;
extern double dz_const;

extern double Lx, Ly, depth;

extern double H_lito;
extern double escala_viscosidade;

extern PetscInt particles_per_ele;

PetscErrorCode _DMLocatePoints_DMDARegular_IS(DM dm,Vec pos,IS *iscell)
{
	
	PetscInt p,n,bs,npoints,si,sj,sk,milocal,mjlocal,mklocal,mx,my,mz;
	DM dmregular;
	PetscInt *cellidx;
	PetscScalar *coor;
	PetscReal dx,dy,dz;
	PetscErrorCode ierr;
	PetscMPIInt rank;
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	ierr = VecGetLocalSize(pos,&n);CHKERRQ(ierr);
	ierr = VecGetBlockSize(pos,&bs);CHKERRQ(ierr);
	npoints = n/bs;
	
	//printf("npoints: %d\n",npoints);
	
	PetscMalloc1(npoints,&cellidx);
	
	ierr = DMGetApplicationContext(dm,(void**)&dmregular);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dmregular,&si,&sj,&sk,&milocal,&mjlocal,&mklocal);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dmregular,NULL,&mx,&my,&mz,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
	
	
	
	//printf("ola!!!\n\n");
	
	dx = dx_const;
	dy = dy_const;
	dz = dz_const;
	
	ierr = VecGetArray(pos,&coor);CHKERRQ(ierr);
	for (p=0; p<npoints; p++) {
		PetscReal coorx,coory,coorz;
		PetscInt mi,mj,mk;
		
		coorx = coor[3*p];
		coory = coor[3*p+1];
		coorz = coor[3*p+2];
		
		mi = (PetscInt)( (coorx)/dx );
		mj = (PetscInt)( (coory)/dy );
		mk = (PetscInt)( (coorz+depth)/dz );
		
		
		cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;

		if ((mk >= sk) && (mk < sk + mklocal)){
			if ((mj >= sj) && (mj < sj + mjlocal)) {
				if ((mi >= si) && (mi < si + milocal)) {
					cellidx[p] = (mi-si) + (mj-sj) * milocal + (mk-sk) * milocal * mjlocal;
				}
			}
		}
		if (coorx < 0.0) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
		if (coorx >  Lx) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
		if (coory < 0.0) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
		if (coory >  Ly) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
		if (coorz < -depth) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
		if (coorz > 0.0) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
	}
	ierr = VecRestoreArray(pos,&coor);CHKERRQ(ierr);
	
	ierr = ISCreateGeneral(PETSC_COMM_SELF,npoints,cellidx,PETSC_OWN_POINTER,iscell);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

PetscErrorCode DMLocatePoints_DMDARegular(DM dm,Vec pos,DMPointLocationType ltype, PetscSF cellSF)
{
	IS iscell;
	PetscSFNode *cells;
	PetscInt p,bs,npoints,nfound;
	const PetscInt *boxCells;
	PetscErrorCode ierr;
	
	ierr = _DMLocatePoints_DMDARegular_IS(dm,pos,&iscell);CHKERRQ(ierr);
	
	//PetscPrintf(PETSC_COMM_WORLD,"teste swarm\n");
	
	ierr = VecGetLocalSize(pos,&npoints);CHKERRQ(ierr);
	ierr = VecGetBlockSize(pos,&bs);CHKERRQ(ierr);
	npoints = npoints / bs;
	
	ierr = PetscMalloc1(npoints, &cells);CHKERRQ(ierr);
	ierr = ISGetIndices(iscell, &boxCells);CHKERRQ(ierr);
	
	for (p=0; p<npoints; p++) {
		cells[p].rank  = 0;
		cells[p].index = DMLOCATEPOINT_POINT_NOT_FOUND;
		
		cells[p].index = boxCells[p];
	}
	ierr = ISRestoreIndices(iscell, &boxCells);CHKERRQ(ierr);
	ierr = ISDestroy(&iscell);CHKERRQ(ierr);
	nfound = npoints;
	ierr = PetscSFSetGraph(cellSF, npoints, nfound, NULL, PETSC_OWN_POINTER, cells, PETSC_OWN_POINTER);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

PetscErrorCode DMGetNeighbors_DMDARegular(DM dm,PetscInt *nneighbors,const PetscMPIInt **neighbors)
{
	DM dmregular;
	PetscErrorCode ierr;
	
	ierr = DMGetApplicationContext(dm,(void**)&dmregular);CHKERRQ(ierr);
	ierr = DMGetNeighbors(dmregular,nneighbors,neighbors);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode SwarmViewGP(DM dms,const char prefix[])
{
	PetscReal *array;
	PetscInt *iarray;
	PetscReal *geoq_fac;
	PetscReal *rho_fac;
	PetscInt npoints,p,bs;
	FILE *fp;
	char name[PETSC_MAX_PATH_LEN];
	PetscMPIInt rank;
	PetscErrorCode ierr;
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s-rank%d.txt",prefix,rank);
	fp = fopen(name,"w");
	if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open file %s",name);
	ierr = DMSwarmGetLocalSize(dms,&npoints);CHKERRQ(ierr);
	printf("npoints = %d\n",npoints);
	ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"itag",NULL,NULL,(void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"geoq_fac",NULL,NULL,(void**)&geoq_fac);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"rho_fac",NULL,NULL,(void**)&rho_fac);CHKERRQ(ierr);
	for (p=0; p<npoints; p++) {
		if (iarray[p]==0)
			fprintf(fp,"%+1.4e %+1.4e %+1.4e %d %1.4e %1.4e\n",array[3*p],array[3*p+1],array[3*p+2],iarray[p],(double)geoq_fac[p],(double)rho_fac[p]);
	}
	ierr = DMSwarmRestoreField(dms,"itag",NULL,NULL,(void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"geoq_fac",NULL,NULL,(void**)&geoq_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"rho_fac",NULL,NULL,(void**)&rho_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	fclose(fp);
	PetscFunctionReturn(0);
}




PetscErrorCode createSwarm()
{
	PetscErrorCode ierr=0;

	PetscMPIInt rank;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	
	PetscInt bs,nlocal,p,cont;
	
	PetscReal *array;
	PetscInt *iarray;
	PetscReal *rarray;
	PetscReal *rarray_rho;
	
	PetscRandom rand;
	
	
	ierr = DMShellCreate(PETSC_COMM_WORLD,&dmcell);CHKERRQ(ierr);
	ierr = DMSetApplicationContext(dmcell,(void*)da_Veloc);CHKERRQ(ierr);
	dmcell->ops->locatepoints = DMLocatePoints_DMDARegular;
	dmcell->ops->getneighbors = DMGetNeighbors_DMDARegular;
	PetscPrintf(PETSC_COMM_WORLD,"teste swarm externo\n");
	
	/* Create the swarm */
	ierr = DMCreate(PETSC_COMM_WORLD,&dms);CHKERRQ(ierr);
	ierr = DMSetType(dms,DMSWARM);CHKERRQ(ierr);
	ierr = DMSetDimension(dms,3);CHKERRQ(ierr);
	
	ierr = DMSwarmSetType(dms,DMSWARM_PIC);CHKERRQ(ierr);
	ierr = DMSwarmSetCellDM(dms,dmcell);CHKERRQ(ierr);
	
	/* init fields */
	ierr = DMSwarmRegisterPetscDatatypeField(dms,"itag",1,PETSC_INT);CHKERRQ(ierr);
	ierr = DMSwarmRegisterPetscDatatypeField(dms,"geoq_fac",1,PETSC_REAL);CHKERRQ(ierr);
	ierr = DMSwarmRegisterPetscDatatypeField(dms,"rho_fac",1,PETSC_REAL);CHKERRQ(ierr);
	ierr = DMSwarmFinalizeFieldRegister(dms);CHKERRQ(ierr);
	
	{
		PetscInt si,sj,sk,milocal,mjlocal,mklocal;
		PetscReal *LA_coors;
		Vec coors;
		PetscInt cnt;
		
		ierr = DMDAGetCorners(da_Veloc,&si,&sj,&sk,&milocal,&mjlocal,&mklocal);CHKERRQ(ierr);
		
		printf("%d: %d %d %d %d %d %d\n",rank,milocal,mjlocal,mklocal,si,sj,sk);
		
		ierr = DMGetCoordinates(da_Veloc,&coors);CHKERRQ(ierr);
		/*VecView(coors,PETSC_VIEWER_STDOUT_WORLD);*/
		ierr = VecGetArray(coors,&LA_coors);CHKERRQ(ierr);
		
		ierr = DMSwarmSetLocalSizes(dms,milocal*mjlocal*mklocal*(particles_per_ele+1),4);CHKERRQ(ierr);
		ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
		
		printf("%d: %d\n",rank,nlocal);
		
		ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		
		printf("bs = %d\n",bs);
		
		
		/*if (rank==0){
			for (int cont=0;cont<10;cont++)
				printf("array %f\n",array[cont]);
		}*/
		
		cnt = 0;
		ierr = PetscRandomCreate(PETSC_COMM_SELF,&rand);CHKERRQ(ierr);
		ierr = PetscRandomSetType(rand,PETSCRAND48);CHKERRQ(ierr);
		ierr = PetscRandomSetInterval(rand,-1.0,1.0);CHKERRQ(ierr);
		for (p=0; p<nlocal; p++) {
			PetscReal px,py,pz,rx,ry,rz;
			
			for (cont=0;cont<particles_per_ele;cont++){
				
				ierr = PetscRandomGetValueReal(rand,&rx);CHKERRQ(ierr);
				ierr = PetscRandomGetValueReal(rand,&ry);CHKERRQ(ierr);
				ierr = PetscRandomGetValueReal(rand,&rz);CHKERRQ(ierr);
				
				px = LA_coors[3*p+0] + (0.5*rx+0.5)*dx_const;
				py = LA_coors[3*p+1] + (0.5*ry+0.5)*dy_const;
				pz = LA_coors[3*p+2] + (0.5*rz+0.5)*dz_const;
				
				if ((px>=0) && (px<=Lx) && (py>=0) && (py<=Ly) && (pz>=-depth) && (pz<=0)) {
					array[bs*cnt+0] = px;
					array[bs*cnt+1] = py;
					array[bs*cnt+2] = pz;
					cnt++;
				}
			}
			
		}
		
		ierr = PetscRandomDestroy(&rand);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		ierr = VecRestoreArray(coors,&LA_coors);CHKERRQ(ierr);
		ierr = DMSwarmSetLocalSizes(dms,cnt,4);CHKERRQ(ierr);
		
		ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
		
		ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		for (p=0; p<nlocal; p++) {
			//iarray[p] = (PetscInt)rank;
			iarray[p] = p%particles_per_ele;
		}
		
		ierr = DMSwarmRestoreField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
		
		ierr = DMSwarmGetField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"rho_fac",&bs,NULL,(void**)&rarray_rho);CHKERRQ(ierr);
		for (p=0; p<nlocal; p++){
			if (array[p*3+2]>-H_lito){
				rarray[p] = escala_viscosidade;
			}
			else rarray[p] = 1.0;
			
			if (array[p*3+2]<-0.8*depth + 0.02*depth*PetscCosReal(3.14159*array[p*3+0]/Lx)){
				rarray_rho[p]=-1000.0;//-0.1;
			}
			else rarray_rho[p]=0.0;
			
		}
		ierr = DMSwarmRestoreField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"rho_fac",&bs,NULL,(void**)&rarray_rho);CHKERRQ(ierr);
		
		ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		
	}
	
	ierr = DMView(dms,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	ierr = SwarmViewGP(dms,"step_0");CHKERRQ(ierr);
	
	
	MPI_Barrier(PETSC_COMM_WORLD);
	
	PetscFunctionReturn(0);
	
}

