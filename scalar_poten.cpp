#include<iostream>
#include "scalar_poten.h"

void scalar_poten3D(nn n,iDouble ***psi,double ****phi,iDouble ***Ux,iDouble ***Uy,iDouble ***Uz, hh h, Cin c,int dir)
{// function return (*phi) matrix Nx*Ny*Nz

    // dir=1, if current move on the left, and -1, if current move on the right

     //double mas_err=1;
     int count_int=0;
     double poisson_lf,poisson_rt;
     iDouble D2x,D2y,D2z;

     int flag_exit=0;
     double MaxDiscr=-1.;// discription
     double hmin2k=pow(min(h.x,min(h.y,h.z)),2.);

     double ***phi_ev;
     double inhx2=1./ (h.x*h.x);
     double inhy2=1./ (h.y*h.y);
     double inhz2=1./ (h.z*h.z);
     double inksi=1.0 / (c.kappa*c.sigma);
     double inhxjsi=h.x*c.j_tr/c.sigma*dir;


     phi_ev=constrDouble3D();





     while (flag_exit==0)//((mas_err>eps_phi)  && count_int<MaxNiterScal)
//                  calculate operator FAx,FAy for input parameters  at each grid node
        {
         MaxDiscr=-1.;
         //in the internal points of superconductor !!!!
         for (int i = n.sx+1; i <= n.ex-1; i++)
            for (int k = n.sz+1; k <= n.ez-1; k++)
             {
                     if(!((n.mx<=i && i<=n.ex) && (n.mbz<=k && k<=n.muz)))
                         for (int j = n.sy+1; j <= n.ey-1; j++)
                     {
                         poisson_lf = ((*phi)[i + 1][j][k] - 2.0*(*phi)[i][j][k] + (*phi)[i - 1][j][k]) *inhx2 +
                                      ((*phi)[i][j + 1][k] - 2.0*(*phi)[i][j][k] + (*phi)[i][j - 1][k]) *inhy2+
                                      ((*phi)[i][j][k+1] -   2.0*(*phi)[i][j][k] + (*phi)[i][j][k-1]) *inhz2;

                         D2x=(psi[i+1][j][k]*Ux[i][j][k]-2.0*psi[i][j][k]+psi[i-1][j][k]*conj(Ux[i-1][j][k])) *inhx2;
                         D2y=(psi[i][j+1][k]*Uy[i][j][k]-2.0*psi[i][j][k]+psi[i][j-1][k]*conj(Uy[i][j-1][k]))*inhy2;
                         D2z=(psi[i][j][k+1]*Uz[i][j][k]-2.0*psi[i][j][k]+psi[i][j][k-1]*conj(Uz[i][j][k-1]))*inhz2;

                         poisson_rt = (inksi)*imag(conj(psi[i][j][k]) * (D2x+ D2y+ D2z));

                         phi_ev[i][j][k] = (*phi)[i][j][k] + h.tsc*(poisson_lf -  poisson_rt);

                         if (MaxDiscr<abs(poisson_lf -  poisson_rt))
                             MaxDiscr=abs(poisson_lf -  poisson_rt);
                     }
             }

         //        boundary and interface conditions=====================================

//                  for (int i=n.sx;i<=n.ex;i++)
//                      for (int k = n.sz; k <= n.ez; k++)
//                  {
//                      phi_ev[i][n.sy][k]=phi_ev[i][n.sy+1][k];//
//                      phi_ev[i][n.ey][k]=phi_ev[i][n.ey-1][k];

//                  }
//                  for (int j=n.sy;j<=n.ey;j++)
//                      for (int k = n.sz; k <= n.ez; k++)
//                  {
//                      phi_ev[n.sx][j][k]=phi_ev[n.sx+1][j][k]+inhxjsi;
//                      phi_ev[n.ex][j][k]=phi_ev[n.ex-1][j][k]-inhxjsi;

//                  }

//                  for (int i=n.sx;i<=n.ex;i++)
//                      for (int j = n.sy; j <= n.ey; j++)
//                  {
//                      phi_ev[i][j][n.sz]=phi_ev[i][j][n.sz+1];//
//                      phi_ev[i][j][n.ez]=phi_ev[i][j][n.ez-1];

//                  }
//========================================================================================
                  // plane ZOY
                  for (int j=n.sy;j<=n.ey;j++)
                      {
                          for (int k=n.sz;k<=n.mbz;k++)
                          {
                          phi_ev[n.ex][j][k]=phi_ev[n.ex-1][j][k]-inhxjsi;
                          }
                          for (int k=n.muz;k<=n.ez;k++)
                          {
                          phi_ev[n.ex][j][k]=phi_ev[n.ex-1][j][k]+inhxjsi;
                          }
                      }
                  for (int j=n.sy;j<=n.ey;j++)
                      {
                          for (int k=n.sz;k<=n.ez;k++)
                          {
                          phi_ev[n.sx][j][k]=phi_ev[n.sx+1][j][k];
                          }
                          for (int k=n.mbz;k<=n.muz;k++)
                          {
                           phi_ev[n.mx][j][k]=phi_ev[n.mx-1][j][k];
                          }
                      }
                  //plane XOZ


                  for (int i=n.sx;i<n.ex;i++)
                      for (int k=n.sz;k<=n.ez;k++)
                      {
                          if(!((n.mx<i && i<=n.ex) && (n.mbz<k && k<n.muz)))
                          {
                          phi_ev[i][n.sy][k]= phi_ev[i][n.sy+1][k];
                          phi_ev[i][n.ey][k]= phi_ev[i][n.ey-1][k];
                          }
                      }

                  // plane XOY
                  for (int i=n.sx;i<n.ex;i++)
                      for (int j=n.sy;j<=n.ey;j++)
                      {
                      phi_ev[i][j][n.sz]=phi_ev[i][j][n.sz+1];
                      phi_ev[i][j][n.ez]=phi_ev[i][j][n.ez-1];
                      }
                  for (int i=n.mx;i<n.ex;i++)
                      for (int j=n.sy;j<=n.ey;j++)
                      {
                     phi_ev[i][j][n.muz]=phi_ev[i][j][n.muz+1];
                     phi_ev[i][j][n.mbz]=phi_ev[i][j][n.mbz-1];
                      }

   //========================================================================
                  //check exit condition
                  if (MaxDiscr<hmin2k && MaxDiscr>0 )// && MaxAbs3D(*AxInd,Ax_old)<eps_B0)
                      flag_exit=1;

//                  mas_err=fabs(phi_ev[0][0][0]-(*phi)[0][0][0]);
//                  for (int i = n.sx; i <= n.ex; i++)      //!!!
//                      for (int k = n.sz; k <= n.ez; k++)
//                              if(!((n.mx<i && i<=n.ex) && (n.mbz<k && k<n.muz)))
//                              {
//                                  for (int j = n.sy; j <= n.ey; j++)
//                                      if (mas_err<fabs(phi_ev[i][j][k]-(*phi)[i][j][k]))
//                                          mas_err=fabs(phi_ev[i][j][k]-(*phi)[i][j][k]);
//                              }
                  for (int i = n.sx; i <= n.ex; i++)    //!!!
                      for (int k = n.sz; k <= n.ez; k++)
                            {
                            if(!((n.mx<i && i<=n.ex) && (n.mbz<k && k<n.muz)))
                                for (int j = n.sy; j <= n.ey; j++)
                                    (*phi)[i][j][k]=phi_ev[i][j][k];

                            }

              count_int++;
          }

//if (mas_err>eps_phi)
//        cout<<mas_err<< "Error for scalar potential"<<endl;

cout<<"scalar iter="<<count_int<<endl;


     freeDouble3D(phi_ev);

}

