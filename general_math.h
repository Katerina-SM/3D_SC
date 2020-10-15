#ifndef GENERAL_MATH_H
#define GENERAL_MATH_H

#include "def_struct.h"

static __inline__ unsigned long long rdtsc(void);

void sum(iDouble ***A, iDouble **B,double k);

void sum3D(iDouble ****A, iDouble ***B,double koef);

void sum3DReal(double ****A, double ***B,double koef);

int isNull3D(iDouble ***B);

void fill(iDouble ***B,double val);

void fill3D(iDouble ****B,double val);

void fill3DReal(double ****B,double val);

void abs_val(iDouble **B,double ***abs_B);

void save_voltage3D(double var,int count,char fname[],string dir);

 void save_general_inform3D(Cin c,nn n,hh h,Coef acoef,string dir);


double max_mass(double **mas, int count_s, int count_y);

double max_mass3D(double ***mas);

double max_abs_mass3D(double ***mas);


double max_mas_AbsComp(iDouble **mas);

double max_mas_AbsComp3D(iDouble ***mas);


void save_abs_data(iDouble **var,int count,double time,char fname[],string dir);

void save_double_data(double **var,int count,double time,char fname[],string dir);

void save_complex_data(iDouble **var,int count,double time,char fname[],string dir);


iDouble **read_complex(char fname[]);

double **read_double(char fname[]);

iDouble ***read_complex3D(char fname[]);
double  ***read_double3D(char fname[]);


void save_test();


void save_abs_data3D(iDouble ***var,int count,double time,char fname[],string dir);

void save_double_data3D(double ***var,int count,double time,char fname[],string dir);

 void save_complex_data3D(iDouble ***var,int count,double time,char fname[],string dir);




double MaxAbs(iDouble **A,iDouble **B);
double MaxAbs3D(iDouble ***A,iDouble ***B);
double MaxAbs3DReal(double ***A,double ***B);
int MaxAbs3DReal3Arr(double ***A,double ***A_old,double ***B,double ***B_old,double ***C,double ***C_old, double hmin2);
void copyComp3D(iDouble ***toCopy,iDouble ****Copy);

//=======================================2D=========================================================
iDouble **constrComp();
void freeComp(iDouble **var);
double **constrDouble();
void freeDouble(double **var);

//=======================================3D=========================================================
double ***constrDouble3D(); //3D matrix
iDouble ***constrComp3D();
void freeComp3D(iDouble ***var);
void freeDouble3D(double ***var);
//==================================================================================================

iDouble As_value(double s, double y, double z, double B_ind);
iDouble Ay_value(double s, double y, double z, double B_ind, double a);
iDouble Az_value(double s, double y, double z, double B_ind);

 string dlts (double x);
 void RecoverSolComp(iDouble ***X,int count,char fname[],char out_psi[]);

 void RecoverSolReal(double ***X,int count,char fname[],char out_psi[] );

#endif // GENERAL_MATH_H
