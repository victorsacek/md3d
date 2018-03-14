//#include <petscksp.h>
//#include <petscksp.h>
//#include <petscmath.h>
///#include <petscdmda.h>
///#include <petscdmswarm.h>
#include <petscsf.h>
//#include <petscdm.h>
#include <petscksp.h>
#include <petscdmda.h>
//#include <petscdmshell.h>
#include <petscdmswarm.h>
#include <petsc/private/dmimpl.h>
#include <petscmath.h>
#include "petscsys.h"

typedef struct {
	PetscScalar u;
	PetscScalar v;
	PetscScalar w;
	//PetscScalar p;
} Stokes;

//PetscErrorCode SwarmViewGP(DM dms,const char prefix[]);

extern DM dms;

extern DM da_Veloc;

extern DM da_Thermal;

extern long Nx,Ny,Nz;

extern double dx_const;
extern double dy_const;
extern double dz_const;

extern double Lx, Ly, depth;

extern Vec local_V,Veloc_weight;

extern Vec geoq_cont,local_geoq_cont;

extern PetscInt particles_per_ele;
extern PetscInt cont_particles;

extern long V_NE;



PetscErrorCode moveSwarm(PetscReal dt)
{
	PetscErrorCode ierr=0;
	
	Stokes					***VV;
	ierr = VecZeroEntries(local_V);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(da_Veloc,Veloc_weight,INSERT_VALUES,local_V);
	ierr = DMGlobalToLocalEnd(  da_Veloc,Veloc_weight,INSERT_VALUES,local_V);
	
	ierr = DMDAVecGetArray(da_Veloc,local_V,&VV);CHKERRQ(ierr);
	
	
	
	PetscInt nlocal,bs,p;
	
	PetscReal *array;
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	
	ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	
	PetscReal N_x0[V_NE], N_y0[V_NE], N_z0[V_NE],strain[6],E2_invariant;
	
	
	PetscReal *strain_fac;
	ierr = DMSwarmGetField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	
	for (p=0; p<nlocal; p++) {
		PetscReal cx,cy,cz,vx,vy,vz;
		PetscReal rx,ry,rz,rfac;
		PetscInt i,j,k;
		PetscInt ii,jj,kk;
		
		PetscReal kx,ky,kz,ex,ey,ez;
		
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
		
		
		/////////cumulative strain
		
		
		kx = 2*rx-1;
		ky = 2*ry-1;
		kz = 2*rz-1;
		
		
		PetscInt cont = 0;
		for (ez=-1.;ez<=1.;ez+=2.){
			for (ey=-1.;ey<=1.;ey+=2.){
				for (ex=-1.;ex<=1.;ex+=2.){
					//N[cont]=(1+ex*kx)*(1+ey*ky)*(1+ez*kz)/8.0;
					N_x0[cont]=ex*(1+ey*ky)*(1+ez*kz)/4.0/dx_const;
					N_y0[cont]=(1+ex*kx)*ey*(1+ez*kz)/4.0/dy_const;
					N_z0[cont]=(1+ex*kx)*(1+ey*ky)*ez/4.0/dz_const;
					cont++;
				}
			}
		}
		
		cont=0;
		
		for (ii=0;ii<6;ii++) strain[ii]=0.0;
		E2_invariant=0.0;
		
		for (kk=k;kk<=k+1;kk++){
			for (jj=j;jj<=j+1;jj++){
				for (ii=i;ii<=i+1;ii++){
					strain[0] += N_x0[cont]*VV[kk][jj][ii].u;
					strain[1] += N_y0[cont]*VV[kk][jj][ii].v;
					strain[2] += N_z0[cont]*VV[kk][jj][ii].w;
					
					strain[3] += N_y0[cont]*VV[kk][jj][ii].u + N_x0[cont]*VV[kk][jj][ii].v;
					strain[4] += N_z0[cont]*VV[kk][jj][ii].v + N_y0[cont]*VV[kk][jj][ii].w;
					strain[5] += N_x0[cont]*VV[kk][jj][ii].w + N_z0[cont]*VV[kk][jj][ii].u;
					
					
					cont++;
				}
			}
		}
		
		
		
		
		
		E2_invariant = (strain[0]-strain[1])*(strain[0]-strain[1]);
		E2_invariant+= (strain[1]-strain[2])*(strain[1]-strain[2]);
		E2_invariant+= (strain[2]-strain[0])*(strain[2]-strain[0]);
		E2_invariant/=6.0;
		E2_invariant+= strain[3]*strain[3];
		E2_invariant+= strain[4]*strain[4];
		E2_invariant+= strain[5]*strain[5];
		
		
		
		
		//strain_fac[p]+= dt*PetscSqrtReal(E2_invariant);//original!!!!
		strain_fac[p]= PetscSqrtReal(E2_invariant);//!!!! não é o cumulativo! apenas o instantaneo.
		
		
		//dt=0.01;
		//vx = 0;
		//vy = -cz;
		//vz = cy;
		
		array[3*p  ] += dt * vx;
		array[3*p+1] += dt * vy;
		array[3*p+2] += dt * vz;
		
		
		
		
	}
	
	ierr = DMSwarmRestoreField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	 
	ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	
	ierr = DMSwarmMigrate(dms,PETSC_TRUE);CHKERRQ(ierr);
	
	//PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN-1,"step%d",tk);
	//ierr = SwarmViewGP(dms,prefix);CHKERRQ(ierr);
	
	//ierr = SwarmViewGP(dms,"step1");CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Veloc,local_V,&VV);CHKERRQ(ierr);
	
	//exit(1);
	
	PetscFunctionReturn(0);
	
}



PetscErrorCode Swarm_add_remove()
{
	PetscErrorCode ierr=0;
	
	PetscInt *carray;
	
	
	PetscScalar             ***qq_cont;
	
	ierr = VecSet(geoq_cont,0.0);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);
	
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
	
	
	PetscInt nlocal,bs,p;
	
	PetscReal *array;
	PetscInt *iarray;
	PetscReal *rarray;
	PetscReal *rarray_rho;
	PetscReal *rarray_H;
	PetscReal *strain_fac;
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	
	ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"cont",&bs,NULL,(void**)&carray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"rho_fac",&bs,NULL,(void**)&rarray_rho);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"H_fac",&bs,NULL,(void**)&rarray_H);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	
	PetscInt Mx=0,mx=10000,My=0,my=10000,Mz=0,mz=10000;
	PetscInt       sx,sy,sz,mmx,mmy,mmz;
	
	ierr = DMDAGetCorners(da_Thermal,&sx,&sy,&sz,&mmx,&mmy,&mmz);CHKERRQ(ierr);

	PetscReal cx,cy,cz,dx,dy,dz;
	PetscInt i,j,k;
	PetscReal cx_v[10],cy_v[10],cz_v[10];
	
	for (p=0; p<nlocal; p++) {

		
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
		
		qq_cont[k][j][i] += 1.0;
		
		carray[p] = k*Nx*Ny + j*Nx + i;
		
		if (Mx<i) Mx=i;
		if (mx>i) mx=i;
		
		if (My<j) My=j;
		if (my>j) my=j;
		
		if (Mz<k) Mz=k;
		if (mz>k) mz=k;
	
	}
	
	PetscInt max_particles_per_ele=particles_per_ele+particles_per_ele/10+2;
	PetscInt min_particles_per_ele=particles_per_ele-particles_per_ele/10-2;
	
	PetscInt kk,pp;
	PetscInt ppp[particles_per_ele*50],cont_p;
	
	PetscInt p_remove[particles_per_ele*50],cont_p_remove=0;
	PetscInt p_i[particles_per_ele*50];
	
	
	PetscInt cont_p_add=0;
	PetscReal p_add_coor[particles_per_ele*50*3];
	PetscReal p_add_r[particles_per_ele*50];
	PetscReal p_add_r_rho[particles_per_ele*50];
	PetscReal p_add_r_H[particles_per_ele*50];
	PetscInt p_add_i[particles_per_ele*50];
	PetscReal p_add_r_strain[particles_per_ele*50];
	
	
	
	
	PetscReal dist,dist_p;
	PetscInt chosen;
	
	PetscRandom rand;
	
	PetscReal rx,ry,rz,xx,yy,zz;
	
	ierr = PetscRandomCreate(PETSC_COMM_SELF,&rand);CHKERRQ(ierr);
	ierr = PetscRandomSetType(rand,PETSCRAND48);CHKERRQ(ierr);
	ierr = PetscRandomSetInterval(rand,-1.0,1.0);CHKERRQ(ierr);
	
	
	for (k=mz; k<=Mz; k++){
		for (j=my; j<=My; j++){
			for (i=mx; i<=Mx; i++){
				if (qq_cont[k][j][i]>max_particles_per_ele){
					qq_cont[k][j][i] -= 1.0;
					/*cont_p=0;
					cx = 0.0;
					cy = 0.0;
					cz = 0.0;
					kk = k*Nx*Ny + j*Nx + i;
					for (p=0; p<nlocal; p++){
						if (carray[p]==kk){
							ppp[cont_p]=p;
							cont_p++;
							cx += array[3*p];
							cy += array[3*p+1];
							cz += array[3*p+2];
						}
					}
					cx/=cont_p;
					cy/=cont_p;
					cz/=cont_p;
					
					dist = 1.0E30;
					
					for (p=0;p<cont_p;p++){
						dx = cx - array[ppp[p]*3];
						dy = cy - array[ppp[p]*3+1];
						dz = cz - array[ppp[p]*3+2];
						
						if (dx*dx+dy*dy+dz*dz<dist){
							pp = ppp[p];
							dist = dx*dx+dy*dy+dz*dz;
						}
					}
					p_remove[cont_p_remove]=pp;
					cont_p_remove++;
					*/
					
					/*
					cont_p=0;
					kk = k*Nx*Ny + j*Nx + i;
					for (p=0; p<nlocal; p++){
						if (carray[p]==kk){
							ppp[cont_p]=p;
							cont_p++;
						}
					}
					
					dist = 1.0E40;
					for (p=0; p<cont_p; p++){
						xx = array[ppp[p]*3];
						yy = array[ppp[p]*3+1];
						zz = array[ppp[p]*3+2];
						dist_p=0.0;
						
						for (pp=0; pp<cont_p; pp++){
							dx = xx - array[ppp[pp]*3];
							dy = yy - array[ppp[pp]*3+1];
							dz = zz - array[ppp[pp]*3+2];
							dist_p+=dx*dx+dy*dy+dz*dz;
						}
						if (dist_p<dist){
							dist = dist_p;
							chosen = ppp[p];
						}
					}
					 
					
					p_remove[cont_p_remove]=chosen;
					cont_p_remove++;
					*/
					
					cont_p=0;
					kk = k*Nx*Ny + j*Nx + i;
					for (p=0; p<nlocal; p++){
						if (carray[p]==kk){
							ppp[cont_p]=p;
							cont_p++;
						}
					}
					
					dist = 1.0E40;
					for (p=0; p<cont_p; p++){
						xx = array[ppp[p]*3];
						yy = array[ppp[p]*3+1];
						zz = array[ppp[p]*3+2];
						dist_p=1.0E40;
						
						for (pp=0; pp<cont_p; pp++){
							dx = xx - array[ppp[pp]*3];
							dy = yy - array[ppp[pp]*3+1];
							dz = zz - array[ppp[pp]*3+2];
							if (dist_p>dx*dx+dy*dy+dz*dz && p!=pp)
								dist_p=dx*dx+dy*dy+dz*dz;
						}
						if (dist_p<dist){
							dist = dist_p;
							chosen = ppp[p];
						}
					}
					
					p_remove[cont_p_remove]=chosen;
					cont_p_remove++;
					
					
					if (cont_p_remove>particles_per_ele*50){
						printf("MUITO\n");
						exit(1);
					}
					
					
					//printf("REMOVEU %d %d %d: %ld!\n",k,j,i,(long)qq_cont[k][j][i]);
				}
				
				if (qq_cont[k][j][i]<min_particles_per_ele){
					qq_cont[k][j][i] += 1.0;
					cont_p=0;
					
					kk = k*Nx*Ny + j*Nx + i;
					for (p=0; p<nlocal; p++){
						if (carray[p]==kk){
							ppp[cont_p]=p;
							cont_p++;
						}
					}
					for (pp=0;pp<10;pp++){
						ierr = PetscRandomGetValueReal(rand,&rx);CHKERRQ(ierr);
						ierr = PetscRandomGetValueReal(rand,&ry);CHKERRQ(ierr);
						ierr = PetscRandomGetValueReal(rand,&rz);CHKERRQ(ierr);
						
						cx_v[pp] = i*dx_const + (0.5*rx+0.5)*dx_const;
						cy_v[pp] = j*dy_const + (0.5*ry+0.5)*dy_const;
						cz_v[pp] = k*dz_const - depth + (0.5*rz+0.5)*dz_const;
						
					}
					
					dist = 0;
					int p_prox,p_prox_total;
					for (pp=0;pp<10;pp++){
						cx = cx_v[pp];
						cy = cy_v[pp];
						cz = cz_v[pp];
						dist_p = 1.0E30;
						for (p=0;p<cont_p;p++){
							dx = cx - array[ppp[p]*3];
							dy = cy - array[ppp[p]*3+1];
							dz = cz - array[ppp[p]*3+2];
							
							if (dx*dx+dy*dy+dz*dz<dist_p){
								p_prox = ppp[p];
								dist_p = dx*dx+dy*dy+dz*dz;
							}
						}
						if (dist<dist_p){
							p_prox_total = p_prox;
							dist=dist_p;
							chosen = pp;
						}
					}
					p_add_coor[cont_p_add*3] = cx_v[chosen];
					p_add_coor[cont_p_add*3+1] = cy_v[chosen];
					p_add_coor[cont_p_add*3+2] = cz_v[chosen];
					
					p_add_i[cont_p_add] = cont_particles%particles_per_ele;
					cont_particles++;
					
					p_add_r[cont_p_add] = rarray[p_prox_total];
					p_add_r_rho[cont_p_add] = rarray_rho[p_prox_total];
					p_add_r_H[cont_p_add] = rarray_H[p_prox_total];
					p_add_r_strain[cont_p_add] = strain_fac[p_prox_total];
					
					//printf("ADDED %d %d %d: !\n",k,j,i);
					//printf("ADDED %lf %lf %lf: !\n",cx_v[chosen],cy_v[chosen],cz_v[chosen]);
					
					
					cont_p_add++;
					if (cont_p_add>particles_per_ele*50){
						printf("MUITO\n");
						exit(1);
					}
				}
				 
			}
		}
	}
	
	
	ierr = PetscRandomDestroy(&rand);CHKERRQ(ierr);
	
	//printf("%d %d   %d %d   %d %d\n",mx,Mx,my,My,mz,Mz);
	
	//printf("b: %d %d   %d %d   %d %d\n",sx,sx+mmx-1,sy,sy+mmy-1,sz,sz+mmz-1);
	
	
	
	ierr = DMSwarmRestoreField(dms,"cont",&bs,NULL,(void**)&carray);CHKERRQ(ierr);
	
	ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"rho_fac",&bs,NULL,(void**)&rarray_rho);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"H_fac",&bs,NULL,(void**)&rarray_H);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	printf("nlocal_%d_antes\n",nlocal);
	
	
	for (pp=0; pp<cont_p_remove; pp++){
		p_i[pp]=pp;
	}
	
	if (cont_p_remove>0)
		PetscSortIntWithPermutation(cont_p_remove,p_remove,p_i);
	
	//for (pp=0; pp<cont_p_remove; pp++){
	for (pp=cont_p_remove-1; pp>=0; pp--){
		DMSwarmRemovePointAtIndex(dms,p_remove[p_i[pp]]);
	}
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	printf("nlocal_%d %d %d_depois\n",nlocal,cont_p_remove,particles_per_ele*50);
	
	if (cont_p_add>0){
		ierr = DMSwarmAddNPoints(dms,cont_p_add);
		
		ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"rho_fac",&bs,NULL,(void**)&rarray_rho);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"H_fac",&bs,NULL,(void**)&rarray_H);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
		
		for (pp=0; pp<cont_p_add; pp++){
			array[(nlocal+pp)*3] = p_add_coor[pp*3];
			array[(nlocal+pp)*3+1] = p_add_coor[pp*3+1];
			array[(nlocal+pp)*3+2] = p_add_coor[pp*3+2];
			
			rarray[nlocal+pp] = p_add_r[pp];
			
			rarray_rho[nlocal+pp] = p_add_r_rho[pp];
			
			rarray_H[nlocal+pp] = p_add_r_H[pp];
			
			strain_fac[nlocal+pp] = p_add_r_strain[pp];
			
			iarray[nlocal+pp] = p_add_i[pp];
		}
		
		ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"rho_fac",&bs,NULL,(void**)&rarray_rho);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"H_fac",&bs,NULL,(void**)&rarray_H);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
		
	}
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	printf("nlocal_%d %d %d_depois2\n",nlocal,cont_p_add,particles_per_ele*50);
	
	
	//exit(1);
	PetscFunctionReturn(0);
}

