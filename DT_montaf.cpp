#include <petscksp.h>

//double calc_visco_real(long t);
//void			montaKeVeloc_simplif(PetscReal *Ke, long t,PetscReal *KeG, Vec visc_vec,PetscInt *Hexa_thermal);
PetscErrorCode  montaKeVeloc_simplif(PetscReal *Ke,PetscReal *KeG,PetscReal *Temper_ele, PetscReal *geoq_ele);

void montafeVeloc(PetscReal *fMe);
void montaCeVeloc(PetscReal *Ce);

extern long V_GT;
extern long V_GN;
extern long V_NE;

extern long n_quad_top;

extern Vec fh_DT;

extern PetscReal *Vfe_DT;
extern PetscInt *Indices_Vfe_DT;

extern PetscReal *Vfe;

extern PetscReal *VfMe;
extern PetscReal *Ke_veloc;

extern PetscReal *Ke_veloc_general;

extern PetscReal *VCe;

extern Vec v_vec;
extern Vec p_vec;
extern Vec drho_vec;

extern PetscInt *Hexa_thermal;


extern double dx_const;
extern double dy_const;
extern double dz_const;

extern Vec visc_vec;


PetscErrorCode DT_montaf()
{
	PetscErrorCode ierr=0;
	
	long t,i,j,c2;
	
	PetscReal val_aux;
	PetscInt ind_prov;
	
	ierr = VecSet(fh_DT,0.0);CHKERRQ(ierr);
	
	montafeVeloc(VfMe);
	
	montaCeVeloc(VCe,dx_const,dy_const,dz_const);
	
	double volume = dx_const*dy_const*dz_const;
	//double visco_real;
	
	for (t=0; t<n_quad_top; t++){
		
		//visco_real = calc_visco_real(t);
		
		//montaKeVeloc_simplif(Ke_veloc,t,Ke_veloc_general,visc_vec,Hexa_thermal);
		ierr = montaKeVeloc_simplif(Ke_veloc,Ke_veloc_general,visc_vec,Hexa_thermal);
		
		for (i=0;i<4;i++) Vfe_DT[i]=0.0;
		
		for (i=0;i<4;i++){
			for (j=0;j<V_NE;j++){
				for (c2=0;c2<V_GN;c2++){
					ind_prov = Hexa_thermal[t*V_NE+j]*V_GN+c2;
					VecGetValues(v_vec,1,&ind_prov,&val_aux);
					Vfe_DT[i]+=volume*(Ke_veloc[((i+4)*V_GN+2)*V_GT+j*V_GN+c2]) *
										val_aux;
				}
			}
		}
		
		
		for (i=0;i<V_GT;i++){
			Vfe[i]=0.0;
		}
		for (i=0;i<V_GT;i++){
			for (j=0;j<V_NE;j++){
				ind_prov = Hexa_thermal[t*V_NE+j];
				VecGetValues(drho_vec,1,&ind_prov,&val_aux);
				Vfe[i]+=VfMe[i*V_NE+j]*val_aux;
			}
		}
		for (i=0;i<V_GT;i++) Vfe[i]*=volume;
		
		
		for (i=0;i<4;i++) Vfe_DT[i]-=Vfe[(i+4)*V_GN+2];
		
		
		ind_prov = t;
		VecGetValues(p_vec,1,&ind_prov,&val_aux);
		for (i=0;i<4;i++) {
			Vfe_DT[i]+=-volume*VCe[(i+4)*V_GN+2]*val_aux;
		}
		
		
		for (i=0;i<4;i++) Indices_Vfe_DT[i]=Hexa_thermal[t*V_NE+i+4];
		//fh_DN[Hexa_thermal[t][i+4]]-=Vfe[(i+4)*V_GN+2];
		
		
		VecSetValues(fh_DT,4,Indices_Vfe_DT,Vfe_DT,ADD_VALUES);
		
	}
	
	
	
	return (ierr);
}