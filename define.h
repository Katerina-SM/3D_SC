
#ifndef DEFINE_H
#define DEFINE_H

#include<complex>

#define Nx 10
#define Ny 3
#define Nt 500000
#define Niter 10000



using namespace std;
typedef complex<double> iDouble;

const iDouble imI(0,1);

const double ax=16;// size specimen
const double ay=16;// in units of lambda (2*xi - vortex size )

const double eps_psi=1e-5;
const double eps_B0=2e-6;


#endif // DEFINE_H
