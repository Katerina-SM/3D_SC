#ifndef GLEQ_H
#define GLEQ_H


#include "general_math.h"

void sc_current(double ****Jx,double ****Jy,double ****Jz, iDouble ***psi,iDouble ***Ux,iDouble ***Uy, iDouble ***Uz,Cin c, nn n, hh h);

void GLeq3D(  nn n,iDouble ****psi,double ***phi,iDouble ***Ux,iDouble ***Uy,iDouble ***Uz, hh h, Cin c,double *dis_psi);

#endif // GLEQ_H
