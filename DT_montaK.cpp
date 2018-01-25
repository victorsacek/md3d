#include <petscksp.h>

void montahKe(PetscReal *Ke_DT);

extern PetscReal *hKe;

extern PetscInt *Indices_Kh_DT_I;
extern PetscInt *Indices_Kh_DT_J;

extern Mat Kh_DT;

extern double dx_const;
extern double dy_const;

extern double r06;
extern double r8p9;
extern double r5p9;

extern long n_quad_top;
extern long nodes_top;

extern PetscInt *Quad_mesh;

extern KSP DT_ksp;

extern int verif_Kh_DT;

PetscErrorCode DT_montaK()
{
	PetscErrorCode ierr=0;
	
	long t,ii,pos_aux;
	
	//// 9*1: maximo de ate 9 nos ao redor
	//// de um ponto no interior do modelo vezes
	//// o numero de graus de liberdade por no.
	
	
	if (verif_Kh_DT==0) {
		MatCreateSeqAIJ(PETSC_COMM_SELF,nodes_top,nodes_top,9*1,0,&Kh_DT);
		verif_Kh_DT=1;
	}
	else{
		MatZeroEntries(Kh_DT);
	}
	
	montahKe(hKe);
	
	/*for (ii=0;ii<4;ii++){
		for (jj=0;jj<4;jj++){
			printf("%lf ",hKe[ii*4+jj]);
		}
		printf("\n");
	}*/
	
	
	
	
	
	for (t=0; t<n_quad_top; t++){
		
		//printf("t: %ld %ld",t,n_quad_top);
		
		for (ii=0;ii<4;ii++){
			pos_aux = Quad_mesh[t*4+ii];
			Indices_Kh_DT_I[ii]=pos_aux;
			Indices_Kh_DT_J[ii]=pos_aux;
		}
		
		
		ierr = MatSetValues(Kh_DT,4,Indices_Kh_DT_I,4,Indices_Kh_DT_J,hKe,ADD_VALUES);
		CHKERRQ(ierr);
		
	}
	
	
	printf("pre Assembly\n");
	ierr = MatAssemblyBegin(Kh_DT,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Kh_DT,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	printf("post Assembly\n");
	
	
	//ierr = MatView(Kh_DT,PETSC_VIEWER_STDOUT_SELF);
	
	
	
	KSPSetOperators(DT_ksp,Kh_DT,Kh_DT);
	KSPSetInitialGuessNonzero(DT_ksp,PETSC_TRUE);
	KSPSetFromOptions(DT_ksp);
	
	
	//KSPSetType(DT_ksp,KSPBCGS);
	
	
	
	return (ierr);
}

void montahKe(PetscReal *Ke_DT)
{
	double N[4];
	
	
	double area = dx_const*dy_const;
	
	double kx,ky;
	
	double ex,ey;
	
	long cont;
	
	long i,j;
	
	for (i=0;i<4;i++){
		for (j=0;j<4;j++){
			Ke_DT[i*4+j]=0.0;
		}
	}
	
	double Hx,Hy,prodH;
	
	for (kx=-r06; kx<=r06; kx+=r06){
		if (kx==0) Hx=r8p9;
		else Hx=r5p9;
		for (ky=-r06; ky<=r06; ky+=r06){
			if (ky==0) Hy=r8p9;
			else Hy=r5p9;
			
			prodH = Hx*Hy;
			cont=0;
			for (ex=-1.;ex<=1.;ex+=2.){
				for (ey=-1.;ey<=1.;ey+=2.){
					N[cont]=(1+ex*kx)*(1+ey*ky)/4.0;
					cont++;
				}
			}
			
			
			for (i=0;i<4;i++){
				for (j=0;j<4;j++){
					//o fator de 2 vem do Jacobiano
					//(1/4 para area e 1/8 para volume, ou seja
					//o fator desta integral de superficie eh o dobro da
					//integral de volume)
					Ke_DT[i*4+j]+=2*prodH*N[i]*N[j]*area;
				}
			}
		}
	}
	
	
	
	
}