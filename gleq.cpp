#include "stdio.h"
#include <string.h>
#include "gleq.h"


void sc_current(double ****Jx,double ****Jy,double ****Jz, iDouble ***psi,iDouble ***Ux,iDouble ***Uy, iDouble ***Uz,Cin c, nn n, hh h)
{

    //    calculate J.................................................................................
    double inkx,inky,inkz;
    inkx=1/(c.kappa*h.x);
    inky=1/(c.kappa*h.y);
    inkz=1/(c.kappa*h.z);

    // aproximated solutions with an accuracy O(h^2)




    for(int i=n.sx;i<=n.ex;i++)
        for (int j=n.sy;j<=n.ey;j++)
            for (int k=n.sz;k<=n.ez;k++)
                {
                 if(!((n.mx<i && i<=n.ex) && (n.mbz<k && k<n.muz)))
                    {
                    if (i<=n.ex-1)
                             (*Jx)[i][j][k]=inkx*imag(Ux[i][j][k]*(psi[i+1][j][k])*conj(psi[i][j][k]));
                    if (j<=n.ey-1)
                             (*Jy)[i][j][k]=inky*imag(Uy[i][j][k]*(psi[i][j+1][k])*conj(psi[i][j][k]));
                    if (k<=n.ez-1)
                             (*Jz)[i][j][k]=inkz*imag(Uz[i][j][k]*(psi[i][j][k+1])*conj(psi[i][j][k]));
                    }
                }




}



void GLeq3D(  nn n,iDouble ****psi,double ***phi,iDouble ***Ux,iDouble ***Uy,iDouble ***Uz, hh h, Cin c,double *dis_psi)
 {  // function solve  Ginzburg Landau equation
    // function return Jx,Jy,psi

    iDouble ***FPsi;   //FPsi=zeros(Nx+1,Ny+1)+0i;


   iDouble cyclic;
   double inhx2,inhy2,inhz2;
   iDouble ik;


   inhx2= 1/(pow(c.kappa*h.x,2.0));
   inhy2= 1/(pow(c.kappa*h.y,2.0));
   inhz2= 1/(pow(c.kappa*h.z,2.0));
   ik=imI*c.kappa;



    FPsi=constrComp3D();


    // calculate operator Fpsi  for input parameters psi,U
    // at each grid node inside the superconductor region

// we solve equation in the internal points of superconductor
        for(int i=n.sx+1;i<=n.ex-1;i++)
             for (int k=n.sz+1;k<=n.ez-1;k++)
                 {
                 if(!((n.mx<=i && i<=n.ex) && (n.mbz<=k && k<=n.muz)))
                    for (int j=n.sy+1;j<=n.ey-1;j++)
                        {

                        cyclic=0;

                        cyclic+=inhx2* (-conj(Ux[i-1][j][k])*(*psi)[i-1][j][k]+2.0*(*psi)[i][j][k]-Ux[i][j][k]*(*psi)[i+1][j][k]) ;
                        cyclic+=inhy2* (-conj(Uy[i][j-1][k])*(*psi)[i][j-1][k]+2.0*(*psi)[i][j][k]-Uy[i][j][k]*(*psi)[i][j+1][k]) ;
                        cyclic+=inhz2* (-conj(Uz[i][j][k-1])*(*psi)[i][j][k-1]+2.0*(*psi)[i][j][k]-Uz[i][j][k]*(*psi)[i][j][k+1]) ;



                        FPsi[i][j][k]=(1-pow(abs((*psi)[i][j][k]),2.0))*((*psi)[i][j][k])-
                                        cyclic-ik*(phi[i][j][k])*((*psi)[i][j][k]);


                        }
                }

        for(int i=n.sx+1;i<=n.ex-1;i++)
             for (int k=n.sz+1;k<=n.ez-1;k++)
                 {
                 if(!((n.mx<=i && i<=n.ex) && (n.mbz<=k && k<=n.muz)))
                    for (int j=n.sy+1;j<=n.ey-1;j++)
                     {
                        (*psi)[i][j][k]+=h.t*FPsi[i][j][k];//psi=psi+h.t*FPsi;
                      }
                  }
     //save_complex_data3D(FPsi,0.,0.,"FPsi_");
      *dis_psi=max_mas_AbsComp3D(FPsi);

   // we correct value at the boundary of superconductor
   // boundary conditions for psi.....................................................................
        // plane ZOY
        for (int j=n.sy;j<=n.ey;j++)
            {
                for (int k=n.sz;k<=n.mbz;k++)
                    (*psi)[n.ex][j][k]=0.;

                for (int k=n.muz;k<=n.ez;k++)
                    (*psi)[n.ex][j][k]=0.;

            }
        for (int j=n.sy;j<=n.ey;j++)
            {
                for (int k=n.sz;k<=n.ez;k++)
                    (*psi)[n.sx][j][k]=Ux[n.sx][j][k]*((*psi)[n.sx+1][j][k]);    //(*psi)[n.ex][j][k]=

                for (int k=n.mbz;k<=n.muz;k++)
                    (*psi)[n.mx][j][k]=conj(Ux[n.mx-1][j][k])*((*psi)[n.mx-1][j][k]);

            }

        //plane XOZ
        for (int i=n.sx;i<n.ex;i++)
            for (int k=n.sz;k<=n.ez;k++)
            {
                if(!((n.mx<i && i<=n.ex) && (n.mbz<k && k<n.muz)))
                {
                (*psi)[i][n.sy][k]=Uy[i][n.sy][k]*((*psi)[i][n.sy+1][k]);
                (*psi)[i][n.ey][k]=conj(Uy[i][n.ey-1][k])*(*psi)[i][n.ey-1][k];
                }
            }

        // plane XOY
        for (int i=n.sx;i<n.ex;i++)
            for (int j=n.sy;j<=n.ey;j++)
            {
            (*psi)[i][j][n.sz]=Uz[i][j][n.sz]*((*psi)[i][j][n.sz+1]);
            (*psi)[i][j][n.ez]=conj(Uz[i][j][n.ez-1])*(*psi)[i][j][n.ez-1];
            }
        for (int i=n.mx;i<n.ex;i++)
            for (int j=n.sy;j<=n.ey;j++)
            {
            (*psi)[i][j][n.muz]=Uz[i][j][n.muz]*((*psi)[i][j][n.muz+1]);
            (*psi)[i][j][n.mbz]=conj(Uz[i][j][n.mbz-1])*(*psi)[i][j][n.mbz-1];
            }




        freeComp3D(FPsi);

}

