

long Nx,Ny,Nz;
long layers;

long T_NE = 8;
long T_GN = 1;

long DIMEN = 3;

long GaussQuad = 27;


long V_NE = 8;
long V_GN = 3;
long V_GT = V_NE*V_GN;

double Lx, Ly, depth;

/////////

double seg_per_ano = 365.0*24.0*3600.0;

double dt_calor = 40000.0;
double dt_calor_sec=dt_calor*seg_per_ano;

double tempo;


double alpha_thermal=0.5;
double comp_alpha_thermal = 1.0 - alpha_thermal;





/////////

double dx_const;
double dy_const;
double dz_const;

int ContMult;

long stepMAX;
double timeMAX;
double dt_MAX;

long print_step;

double visco_r;

double visc_MAX;
double visc_MIN;

double escala_viscosidade;

double veloc_superf;

double RHOM;
double alpha_exp_thermo;
double kappa;


double gravity;

double Delta_T;

double H_lito;

double H_per_mass;
double c_heat_capacity;

int T_initial_cond;
int rheol;

double beta_max;
double ramp_begin;
double ramp_end;

int bcv_top_normal;
int bcv_top_slip;

int bcv_bot_normal;
int bcv_bot_slip;

int bcv_left_normal;
int bcv_left_slip;

int bcv_right_normal;
int bcv_right_slip;

int bcT_top;

int bcT_bot;

int bcT_left;

int bcT_right;

///////

PetscReal *TKe, *TCe, *TFe, *TCe_fut, *TMe, *Ttotal, *Ttotal_b;


PetscReal *T_vec_aux_ele;
PetscReal *T_vec_aux_ele_final;

double r06 = 0.7745966692414834; //sqrt(0.6)
double r8p9 = 8.0/9.0;
double r5p9 = 5.0/9.0;


PetscReal *NT;
PetscReal *NT_x;
PetscReal *NT_y;
PetscReal *NT_z;


////////

Vec v_vec;
Vec v_vec_fut;

PetscInt *indice_aux_vec_ele;

PetscReal *v_vec_aux_ele;

/////////

Mat TA, TB;
Vec Tf, Temper;

Vec dRho;

DM da_Thermal;

KSP T_ksp;


/////////

Mat VA, VB, VG;
Vec Vf, Veloc, Veloc_fut;

Vec Pressure;

DM da_Veloc;

KSP V_ksp;

PetscReal *Vfe;

PetscReal *Ke_veloc;
PetscReal *Ke_veloc_final;

PetscReal *Ke_veloc_general;

PetscReal *VCe;
PetscReal *VfMe;

PetscInt Verif_VG=0;


Vec rk_vec2;

Vec rk_vec;
Vec sk_vec;
Vec gs_vec;
Vec uk_vec;

Vec Veloc_Cond;



