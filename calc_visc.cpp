#include <petscksp.h>
//#include <stdio.h>
#include <math.h>

extern int rheol;

extern double visc_MAX;
extern double visc_MIN;

extern int geoq_on;

extern double visco_r;

extern double Delta_T;

double calc_visco_ponto(double T,double z,double geoq_ponto){
	
	double visco_real;
	
	if (rheol==0)	visco_real = visco_r;
	
	if (rheol==1){
		double r=20.0;
		double Q = 225.0/log(r)-0.25*log(r);
		double G = 15./log(r)-0.5;
		
		return(geoq_ponto*visco_r*exp(  Q/(T/Delta_T+G) - Q/(0.5+G)    ));
	}
	
	if (rheol==2){
		double R = 8.31;     // J/mol/K
		double E = 120000.0; // J/mol
		
		
		double Tb = Delta_T+273.0;
		
		double aux = E*(1.0/(T+273.0)-1.0/Tb)/R;
		visco_real = visco_r*exp(aux);
		
	}
	
	
	if (rheol==3){
		double R = 8.31;     // J/mol/K
		double E = 120000.0; // J/mol
		
		
		double Tb = 1300.0+273.0;
		
		double aux = E*(1.0/(T+273.0)-1.0/Tb)/R;
		visco_real = visco_r*exp(aux);
		
	}
	
	if (rheol==4){
		double R = 8.31;     // J/mol/K
		double E = 120000.0; // J/mol
		
		
		double Tb = Delta_T+273.0;
		
		double aux = E*(1.0/(T+273.0)-1.0/Tb)/R;
		visco_real = visco_r*exp(aux);
		
	}
	
	if (rheol==5){
		double R = 8.31;     // J/mol/K
		double E = 240000.0; // J/mol
		
		double b = 1.0E7;
		
		
		double Tb = Delta_T+273.0;
		
		double aux = -(T+273)*E/(R*Tb*Tb);
		visco_real = visco_r*b*exp(aux);
		
	}
	
	if (rheol==6){
		double R = 8.31;     // J/mol/K
		double E = 240000.0; // J/mol
		
		double b = 1.0E7;
		
		
		double Tb = Delta_T+273.0;
		
		double aux = -(T+273)*E/(R*Tb*Tb);
		visco_real = visco_r*b*exp(aux);
		
		
	}
	
	
	if (rheol>6){
		printf("rheol error: larger than maximum available option\n");
		exit(1);
	}
	
	if (geoq_on)
		visco_real *= geoq_ponto;
	
	
	if (visco_real>visc_MAX) visco_real=visc_MAX;
	if (visco_real<visc_MIN) visco_real=visc_MIN;
	
	return(visco_real);
	
	return(0);
	
	
	
}
