#ifndef DEF_STRUCT_H
#define DEF_STRUCT_H

#include<complex>
#include<vector>

#define Nx 57//87                   //240
#define Ny 57//111                   //350
#define Nz 25//107//7

#define Nt 990000 // General iteration number
#define Niter 5000//iteration number  due to find stable GL solution

#define NiterScal 200000     // maximal iteration number to solve GL equation at the first step before solving equation for scalar potential
#define NiterVect 10000     // maximal iteration number to solve GL equation at the first step before solving equation for vector potential

#define MaxNiterVect 200000// maximal iteration number to solve vector potential equation
#define MaxNiterScal 200000// maximal iteration number to solve scalar potential equation
#define PI 3.141592653589




using namespace std;
typedef complex<double> iDouble;

const iDouble imI(0,1);
const iDouble im0(0, 0);

const double ax=3.012688;//3.4059;// size of general domain where the spieces are situated
const double ay=3.012688;;//0.6;//4.436; // in units of lambda (2*xi - vortex size )
const double az=1.321344;//4.198; // thicknes of domain

const double eps_psi=1e-4;
const double eps_B0=1e-5;
const double eps_phi=2e-5;
const double eps_discr=2e-5;





const double l_mfp=6.0*pow(10,-9);			//5.7nm for NanoLetters
const double c = 3.0 * pow(10, 8);		//light velocity in m/s
const double vF = 6.0*pow(10, 5);			//Fermi velocity in m/s
const double D = l_mfp*vF/3.0;			//diffusion coefficient in m^2/s
const double kT = 0.95;					//dimensionless, equal T/Tc
//const double nV = 5.6*pow(10, 28);		//valency electrons concentration in 1/m^3
//const double C_e = 0.281*pow(10, -7);		//q^2_e/m_e in Si (Couloumb/kg)
const double C_fi = 0.11*pow(10, -17);	//equal Plank constant divide 2x electron charge [erg/SGS]
//const double C_force = 9.08*pow(10, 20);	//dimensionless, equal (c*Hc*m_s)/(2*sqrt(2)*PI*e_s*h_planck) - coeff from LorForce
//const double sigma_dim = l_mfp/(3.72*pow(10, -16)); //conductivity in 1/(Om*m)	//Drude model: C_e*pow(kT, 4)*l_mfp*nV / vF;
const double lambda_0 = 39.0*pow(10, -9);	//const for lambda
const double xi_0 = 39.0*pow(10, -9);		//const for xi
const double F_0 = 2.068*pow(10, -11);	//magnetic flow quant in Gauss*m^2
//const double psi_0_square = 1.81*pow(10,20);								//psi^2_0 from equations
const double lambda = lambda_0*sqrt(xi_0 / (2.0*(1.0 - kT)*1.33*l_mfp));	//London's penetration depth in m
const double xi = 0.855*sqrt(xi_0*l_mfp / (1.0 - kT));					//coherency length in m
const double kappa = lambda / xi;											//GL parameter
const double tau = pow(xi, 2) / D;										//characteristic time in s
const double Hc = F_0 / (2.0*PI*lambda*xi*sqrt(2.0));						//magnetic quant unit in Gauss
//const double B_0 = sqrt(2.0)*Hc*pow(10, -4);								//magnetic field unit in Tesla
const double fi_0 = pow(10,6)*300.0*kappa*D*C_fi/pow(xi,2);				//voltage unit in muV
const double j_0 = (3.34*pow(10,-6))*F_0*c/(8.0*PI*PI*pow(lambda,2)*xi);	//current density unit
const double sigma_0 = pow(c, 2)/(4.0*PI*kappa*kappa*D*8.988*pow(10,9));	//conductivity unit in 1/(Om*m)
//const double delta = 60.0*pow(10, -9);	//cut width in m		(59nm for NanoLetters)
//const double L = 5.*pow(10, -6);			//cylinder length in m	(3.5nm for NanoLetters)
//const double R_nondim = R/lambda;         //dimensionless radius
//const double delta_nondim = delta/lambda; //dimensionless cut
//const double a = (2.*PI*R - delta)/lambda;//dimensionless S size	(a = 6.097 - length)
//const double b = L / lambda;				//dimensionless Y size	(b = 12.544 - width)
//const double sigma = sigma_dim / sigma_0;		//conductivity dimensionless
//const double hs = a/Ns;
//const double hy = b/Ny;
//const double h = max(hs,hy);
//const double ht = 1.0*h*h;
//const double h_tau = ht/10.;
//const double  eps = 0.00002;
//const double  N= max(Ns, Ny);

//ostringstream directoryCurr;


struct hh
{
    double x;
    double y;
    double z;
    double t;
    double tau;
    double tsc;
 };

struct nn
{
    int sx;
    int mx;
    int ex;

    int sy;
    int ey;

    int sz;
    int mbz;
    int muz;
    int ez;
};
struct Cin
{
    double kappa;
    double B0x;
    double B0y;
    double B0z;
    double sigma;
    double j_tr;

};

struct Coef
{   double a;
    double sx;
    double sy;
 };

struct koefGL
{
    double koef;
    double koef_min;
    double koef_max;
 };



#endif // DEF_STRUCT_H
