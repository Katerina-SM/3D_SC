#include<iostream>
#include "gl_vectorpoten.h"
#include "gleq.h"
#include "scalar_poten.h"
#include "general_math.h"

#include <omp.h>
//#include <cstring>

//#include <sys/types.h>
#include <sys/stat.h>






void Initialise(nn *n, hh  *h, Cin *c, Coef *acoef, koefGL *kGL)
{
    double hmin2;
    double thick=1.879;// thick of plane( in lambda units )

    int hole=3;//  number of level
    int shift=4;
    // order parameter are in nodes (i,j), i=n.sx .. n.ex, j=n.sx .. n.ex

    h->x=ax/Nx;     h->y=ay/Ny;     h->z=az/Nz;

    n->sx=shift;//int(floor(Nx*(0.2))+1);
    n->mx=shift+7;//floor(thick/h->x);
    n->ex=Nx-shift;//int(floor(Nx*(0.8)));



    n->sy=shift;//int(floor(Ny*(0.2)) +1);
    n->ey=Ny-shift;//int(floor(Ny*(0.8)));

    n->sz=2;
    n->mbz=n->sz+7;//floor(thick/h->z);
    n->muz=n->mbz+7;//hole;
    n->ez=n->muz+7;//floor(thick/h->z);//Nz-2;

    // GL parameter
    c->kappa=lambda/xi;


    // magnetic field along OZ
    c->B0z=0.25;//0.25;//0.38;//0.48452237;    // (in units of sqrt(2)*Hc, Hc - magnetic quant unit, Gauss)
    c->B0x=0.;
    c->B0y=0.;

    c->sigma=0.527756; // Normal conductivity

    c->j_tr=0;//0.1;//0.2//0.2816439;//0.3313;


    double hmin=min(h->x,min(h->y,h->z));

    kGL->koef_max=9./(10.*hmin);
    kGL->koef_min=kGL->koef_max/10.;
    kGL->koef=kGL->koef_min;

    hmin2=pow(min(min(h->x,h->y),h->z),2.0);

    h->t = hmin2*kGL->koef;// time for  GL equation
    h->tau=hmin2/6.;// time for  Vector potential equation
    h->tsc=hmin2/6.;//10.;// time for  Scalar potential equation

    acoef->a=1.;
    acoef->sx=Nx*h->x/2;
    acoef->sy=Ny*h->y/2;





}

void CalculateLength(double ****len,hh h)
{
    double koef=1/(4*PI)*h.x*h.y*h.z;
    double lenOp;
    double val;
    double hx2=pow(h.x,2.),hy2=pow(h.y,2.),hz2=pow(h.z,2.);

    for (int i=0;i<=Nx;i++)
        for (int j=0;j<=Ny;j++)
            for (int k=0;k<=Nz;k++)
                {
                val=koef/sqrt(i*(i*hx2)+j*(j*hy2)+k*(k*hz2));

                lenOp=(i==0 && j==0 && k==0)? 0:val;

                (*len)[i][j][k]=lenOp;
            }

}



void SCDomainValOMP(int ii,int jj,int kk,int flagX,int flagY,int flagZ,
                 double ***len,
                 double ***Jx,double  ***Jy,double ***Jz,double *valX,double *valY,double *valZ, nn n)
{
    //double len;
    //double koef=1/(4*PI)*h.x*h.y*h.z;
    double vX=0,vY=0,vZ=0;



    int il,jl,kl;

    *valX=0;*valY=0;*valZ=0;

    // double tst = omp_get_wtime();
     //double t;

//#pragma omp parallel
    {

    (vX)=0;
    (vY)=0;
    (vZ)=0;

   #pragma omp parallel for reduction(+:vX) reduction(+:vY) reduction(+:vZ)
    for(int i=n.sx;i<=n.ex;i++)
        for(int j=n.sy;j<=n.ey;j++)
            for(int k=n.sz;k<=n.ez;k++)
            {
               // We suggest that all points (i,j,k) belong superconductor
               //len=koef/sqrt(pow(double(ii-i)*h.x,2.)+pow(double(jj-j)*h.y,2.)+pow(double(kk-k)*h.z,2.));
                il=abs(ii-i);
                jl=abs(jj-j);
                kl=abs(kk-k);
    
                if (flagX==1)
                    vX +=Jx[i][j][k]*(len[il][jl][kl]);//3D
                if (flagY==1)
                    vY +=Jy[i][j][k]*(len[il][jl][kl]);
                if (flagZ==1)
                    vZ +=Jz[i][j][k]*(len[il][jl][kl]);
    
               }


    }

   //  t = omp_get_wtime() - tst;        //printf("Thread %d execution time: %.6f sec.\n",               omp_get_thread_num(), t);
    // cout<< "time is "<< t<<endl;

    *valX=vX;
    *valY=vY;
    *valZ=vZ;


}

void SCDomainVal(int ii,int jj,int kk,int flagX,int flagY,int flagZ,
                 double ***len,
                 double ***Jx,double  ***Jy,double ***Jz,double *valX,double *valY,double *valZ, nn n)
{
    //double len;
    //double koef=1/(4*PI)*h.x*h.y*h.z;

    int il,jl,kl;

    double vX=0,vY=0,vZ=0;


   // double tst = omp_get_wtime();
   // double t;

for(int i=n.sx;i<=n.ex;i++)
    for(int k=n.sz;k<=n.ez;k++)
            if(!((n.mx<i && i<=n.ex) && (n.mbz<k && k<n.muz)))
                for(int j=n.sy;j<=n.ey;j++)
                    {
                // We suggest that all points (i,j,k) belong superconductor

                    //len=koef/sqrt(pow(double(ii-i)*h.x,2.)+pow(double(jj-j)*h.y,2.)+pow(double(kk-k)*h.z,2.));

                    il=abs(ii-i);
                    jl=abs(jj-j);
                    kl=abs(kk-k);

                    if (flagX==1)
                        vX+=Jx[i][j][k]*(len[il][jl][kl]);//3D
                    if (flagY==1)
                        vY+=Jy[i][j][k]*(len[il][jl][kl]);
                    if (flagZ==1)
                        vZ+=Jz[i][j][k]*(len[il][jl][kl]);
                    }
*valX=vX;
*valY=vY;
*valZ=vZ;

//t = omp_get_wtime() - tst;        //printf("Thread %d execution time: %.6f sec.\n",               omp_get_thread_num(), t);
//cout<< "time  no OMP is "<< t<<endl;
}



void bilinearInterpXOY(double ****Q,int nf_h,int nf_e,int sizeX,int sizeY,int kx,int ky)
{// we suppose that we know  matrix elements Nx*Ny*Nz     x
    double koef=1./(kx*ky);
    int nf;

// down and up plane
    //horisontal
    for (int j=0;j<=sizeY;j+=ky)
        for (int i=0;i<sizeX;i+=kx)
            for( int ii=1;ii<=kx-1;ii++)
                {
                nf=nf_h;
                 (*Q)[i+ii][j][nf]=koef*((*Q)[i][j][nf]*double((kx-ii)*(ky))+(*Q)[i+kx][j][nf]*double(ii*ky));
                nf=nf_e;
                (*Q)[i+ii][j][nf]=koef*((*Q)[i][j][nf]*double((kx-ii)*(ky))+(*Q)[i+kx][j][nf]*double(ii*ky));
                }


    //vertical
    for (int i=0;i<=sizeX;i+=kx)
        for (int j=0;j<sizeY;j+=ky)
               for( int jj=1;jj<=ky-1;jj++)
               {
                   nf=nf_h;
                   (*Q)[i][j+jj][nf]=koef*((*Q)[i][j][nf]*double (kx*(ky-jj))+(*Q)[i][j+ky][nf]*double(kx*jj));
                   nf=nf_e;
                   (*Q)[i][j+jj][nf]=koef*((*Q)[i][j][nf]*double (kx*(ky-jj))+(*Q)[i][j+ky][nf]*double(kx*jj));
               }

    //internal point
    for (int i=0;i<sizeX;i+=kx)
        for (int j=0;j<sizeY;j+=ky)
            for( int ii=1;ii<=kx-1;ii++)
                for( int jj=1;jj<=ky-1;jj++)
                     {
                    nf=nf_h;
                    (*Q)[i+ii][j+jj][nf]=koef*(((*Q)[i][j][nf])*double((kx-ii)*(ky-jj))+(*Q)[i+kx][j][nf]*double((ii)*(ky-jj))+
                        (*Q)[i][j+ky][nf]*double((kx-ii)*(jj))+(*Q)[i+kx][j+ky][nf]*double((ii)*(jj)));
                    nf=nf_e;
                    (*Q)[i+ii][j+jj][nf]=koef*(((*Q)[i][j][nf])*double((kx-ii)*(ky-jj))+(*Q)[i+kx][j][nf]*double((ii)*(ky-jj))+
                        (*Q)[i][j+ky][nf]*double((kx-ii)*(jj))+(*Q)[i+kx][j+ky][nf]*double((ii)*(jj)));
                 }
 }


void bilinearInterpXOZ(double ****Q,int nf_h,int nf_e,int sizeX,int sizeY,int kx,int ky)
{// we suppose that we know  matrix elements Nx*Ny*Nz     x
    double koef=1./(kx*ky);
    int nf;

// down and up plane
    // horisontal(i - correspond i; j - correspond k )
    for (int j=0;j<=sizeY;j+=ky)
        for (int i=0;i<sizeX;i+=kx)
            for( int ii=1;ii<=kx-1;ii++)
                {
                nf=nf_h;
                (*Q)[i+ii][nf][j]=koef*((*Q)[i][nf][j]*double((kx-ii)*(ky))+(*Q)[i+kx][nf][j]*double(ii*ky));
                nf=nf_e;
                (*Q)[i+ii][nf][j]=koef*((*Q)[i][nf][j]*double((kx-ii)*(ky))+(*Q)[i+kx][nf][j]*double(ii*ky));
                }
    //vertical
    for (int i=0;i<=sizeX;i+=kx)
        for (int j=0;j<sizeY;j+=ky)
               for( int jj=1;jj<=ky-1;jj++)
               {nf=nf_h;
                (*Q)[i][nf][j+jj]=koef*((*Q)[i][nf][j]*double (kx*(ky-jj))+(*Q)[i][nf][j+ky]*double(kx*jj));
                nf=nf_e;
                (*Q)[i][nf][j+jj]=koef*((*Q)[i][nf][j]*double (kx*(ky-jj))+(*Q)[i][nf][j+ky]*double(kx*jj));
               }

    //internal point
    for (int i=0;i<sizeX;i+=kx)
        for (int j=0;j<sizeY;j+=ky)
            for( int ii=1;ii<=kx-1;ii++)
                for( int jj=1;jj<=ky-1;jj++)
                     {
                    nf=nf_h;
                    (*Q)[i+ii][nf][j+jj]=koef*(((*Q)[i][nf][j])*double((kx-ii)*(ky-jj))+(*Q)[i+kx][nf][j]*double((ii)*(ky-jj))+
                        (*Q)[i][nf][j+ky]*double((kx-ii)*(jj))+(*Q)[i+kx][nf][j+ky]*double((ii)*(jj)));
                    nf=nf_e;
                    (*Q)[i+ii][nf][j+jj]=koef*(((*Q)[i][nf][j])*double((kx-ii)*(ky-jj))+(*Q)[i+kx][nf][j]*double((ii)*(ky-jj))+
                        (*Q)[i][nf][j+ky]*double((kx-ii)*(jj))+(*Q)[i+kx][nf][j+ky]*double((ii)*(jj)));
                 }
 }

void bilinearInterpYOZ(double ****Q,int nf_h,int nf_e,int sizeX,int sizeY,int kx,int ky)
{// we suppose that we know  matrix elements Nx*Ny*Nz     x
    double koef=1./(kx*ky);
    int nf;

// down and up plane
    // horisontal(i - correspond y; j - correspond z )
     for (int j=0;j<=sizeY;j+=ky)
         for (int i=0;i<sizeX;i+=kx)
            for( int ii=1;ii<=kx-1;ii++)
                {nf=nf_h;
                (*Q)[nf][i+ii][j]=koef*((*Q)[nf][i][j]*double((kx-ii)*(ky))+(*Q)[nf][i+kx][j]*double(ii*ky));
                nf=nf_e;
                (*Q)[nf][i+ii][j]=koef*((*Q)[nf][i][j]*double((kx-ii)*(ky))+(*Q)[nf][i+kx][j]*double(ii*ky));
                }
    //vertical
    for (int i=0;i<=sizeX;i+=kx)
        for (int j=0;j<sizeY;j+=ky)
               for( int jj=1;jj<=ky-1;jj++)
               {nf=nf_h;
                (*Q)[nf][i][j+jj]=koef*((*Q)[nf][i][j]*double (kx*(ky-jj))+(*Q)[nf][i][j+ky]*double(kx*jj));
                nf=nf_e;
                (*Q)[nf][i][j+jj]=koef*((*Q)[nf][i][j]*double (kx*(ky-jj))+(*Q)[nf][i][j+ky]*double(kx*jj));
               }

    //internal point
    for (int i=0;i<sizeX;i+=kx)
        for (int j=0;j<sizeY;j+=ky)
            for( int ii=1;ii<=kx-1;ii++)
                for( int jj=1;jj<=ky-1;jj++)
                     {
                    nf=nf_h;
                    (*Q)[nf][i+ii][j+jj]=koef*(((*Q)[nf][i][j])*double((kx-ii)*(ky-jj))+(*Q)[nf][i+kx][j]*double((ii)*(ky-jj))+
                        (*Q)[nf][i][j+ky]*double((kx-ii)*(jj))+(*Q)[nf][i+kx][j+ky]*double((ii)*(jj)));
                    nf=nf_e;
                    (*Q)[nf][i+ii][j+jj]=koef*(((*Q)[nf][i][j])*double((kx-ii)*(ky-jj))+(*Q)[nf][i+kx][j]*double((ii)*(ky-jj))+
                        (*Q)[nf][i][j+ky]*double((kx-ii)*(jj))+(*Q)[nf][i+kx][j+ky]*double((ii)*(jj)));
                 }
 }

void IntegVectPotenInSpace( double ****AxInt,double ****AyInt,double ****AzInt,double ***len,double ***Jx,double ***Jy,double ***Jz,nn n)
{

//     % numerically we find vector potential which provided by integral formula

//     % input parameters:

//     % (Jx,Jy,Jz) - vector field of current
//     % n - describe boundary of current  and has structure  along OX n.sx,n.ex and along OY n.sy,n.ey
//     % h - grid size (h.x,h.y)
//     % c - unused parameter

//     % var - variant of prepared initial data(in file Initialise)
//     % dimension =2, if we solve in 2D, 3, if we solve in 3D
//     % output parameters:
//     % (Ax,Ay,Az) - vector potential
//     % (Bdx,Bdy,Bdz) - magnetic field  which is provided by differense scheme  rot A, A=(Ax,Ay,Az)
//     % (BxSL,BySL,BzSL) - magnetic field from Bio-Savar law


    double valX=0,valY=0,valZ=0;




    // we determine vector potential  at the boundary from the Integral ratio

//    double search_time;
//    double start_time=clock();


//double t = omp_get_wtime();


              // plane parral XOY
               for (int j=0;j<=Ny;j++)
                for (int k=0;k<=Nz;k++)
                   {
                     SCDomainVal(0,j,k,1,1,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);

                    (*AxInt)[0][j][k]=valX;
                    (*AyInt)[0][j][k]=valY;
                    (*AzInt)[0][j][k]=valZ;

                    SCDomainVal(Nx-1,j,k,1,0,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
                    (*AxInt)[Nx-1][j][k]=valX;

                     SCDomainVal(Nx,j,k,0,1,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
                    (*AyInt)[Nx][j][k]=valY;
                    (*AzInt)[Nx][j][k]=valZ;
                    }
                // plane paralel  XOZ
                   for (int i=0;i<=Nx;i++)
                     for (int k=0;k<=Nz;k++)
                    {
                    SCDomainVal(i,0,k,1,1,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
                    (*AxInt)[i][0][k]=valX;
                    (*AyInt)[i][0][k]=valY;
                    (*AzInt)[i][0][k]=valZ;

                    SCDomainVal(i,Ny,k,1,0,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
                    (*AxInt)[i][Ny][k]=valX;
                    (*AzInt)[i][Ny][k]=valZ;

                    SCDomainVal(i,Ny-1,k,0,1,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
                    (*AyInt)[i][Ny-1][k]=valY;
                    }
             // plane paralel  YOZ
                for (int i=0;i<=Nx;i++)
                 for (int j=0;j<=Ny;j++)

                    {
                     SCDomainVal(i,j,0,1,1,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
                    (*AxInt)[i][j][0]=valX;
                    (*AyInt)[i][j][0]=valY;
                    (*AzInt)[i][j][0]=valZ;

                    SCDomainVal(i,j,Nz,1,1,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
                    (*AxInt)[i][j][Nz]=valX;
                    (*AyInt)[i][j][Nz]=valY;

                    SCDomainVal(i,j,Nz-1,0,0,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
                    (*AzInt)[i][j][Nz-1]=valZ;
                     }

  //t = omp_get_wtime()-t;
  //cout<<"time for all planes "<<t<<endl;

//         unsigned int end_time=clock();
//         search_time=(double(end_time-start_time)/ (double)CLOCKS_PER_SEC);

//         cout<<"calculation of integral boundarz condition takes "<< search_time <<"   sec" <<endl;



}


void IntegVectPotenInSpaceMod( double ****AxInt,double ****AyInt,double ****AzInt,double ***len,double ***Jx,double ***Jy,double ***Jz,nn n)
{

//     % numerically we find vector potential which provided by integral formula

//     % input parameters:

//     % (Jx,Jy,Jz) - vector field of current
//     % n - describe boundary of current  and has structure  along OX n.sx,n.ex and along OY n.sy,n.ey
//     % h - grid size (h.x,h.y)
//     % c - unused parameter

//     % var - variant of prepared initial data(in file Initialise)
//     % dimension =2, if we solve in 2D, 3, if we solve in 3D
//     % output parameters:
//     % (Ax,Ay,Az) - vector potential
//     % (Bdx,Bdy,Bdz) - magnetic field  which is provided by differense scheme  rot A, A=(Ax,Ay,Az)
//     % (BxSL,BySL,BzSL) - magnetic field from Bio-Savar law

// Here we suppose that all components of vector potential have the same size (Nx+1)*(Ny+1)*(Nz+1)

double valX=0,valY=0,valZ=0;
int si=2,sj=2,sk=2;//we suppose that (Nx-1)/si,(Ny-1)/sj,(Nz-1)/sk  - is integer!!!




    // we determine vector potential  at the boundary from the Integral ratio

    double  search_time;
    double  start_time=clock();

    // plane paralel  XOY (top and down planes)
    for (int i=0;i<=Nx;i+=si)
        for (int j=0;j<=Ny;j+=sj)
           {
            SCDomainVal(i,j,0,1,1,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
           (*AxInt)[i][j][0]=valX;
           (*AyInt)[i][j][0]=valY;
           (*AzInt)[i][j][0]=valZ;

           SCDomainVal(i,j,Nz,1,1,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
           (*AxInt)[i][j][Nz]=valX;
           (*AyInt)[i][j][Nz]=valY;

           SCDomainVal(i,j,Nz-1,0,0,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
           (*AzInt)[i][j][Nz-1]=valZ;

            }
//    additional line
    for (int i=0;i<=Nx;i++)
           {
           int j=Ny;
           SCDomainVal(i,j,0,1,0,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
           (*AxInt)[i][j][0]=valX;
           (*AzInt)[i][j][0]=valZ;

           SCDomainVal(i,j,Nz,1,0,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
           (*AxInt)[i][j][Nz]=valX;

           SCDomainVal(i,j,Nz-1,0,0,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
           (*AzInt)[i][j][Nz-1]=valZ;
            }
    for (int j=0;j<=Ny;j++)
           {
            int i=Nx;
            SCDomainVal(i,j,0,0,1,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
           (*AyInt)[i][j][0]=valY;
           (*AzInt)[i][j][0]=valZ;

           SCDomainVal(i,j,Nz,0,1,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
           (*AyInt)[i][j][Nz]=valY;

           SCDomainVal(i,j,Nz-1,0,0,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
           (*AzInt)[i][j][Nz-1]=valZ;
            }

      // plane paralel  ZOY............................................................
       for (int j=0;j<=Ny;j+=sj)
        for (int k=0;k<=Nz-1;k+=sk)//for (int k=0;k<=Nz;k+=sk)
           {
            SCDomainVal(0,j,k,1,1,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);

             (*AxInt)[0][j][k]=valX;
             (*AyInt)[0][j][k]=valY;
             (*AzInt)[0][j][k]=valZ;

            SCDomainVal(Nx,j,k,0,1,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
            (*AyInt)[Nx][j][k]=valY;
            (*AzInt)[Nx][j][k]=valZ;

            SCDomainVal(Nx-1,j,k,1,0,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
            (*AxInt)[Nx-1][j][k]=valX;
            }
       // additional line
       for (int j=0;j<=Ny;j++)
           {int k=Nz;
            SCDomainVal(0,j,k,1,1,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
             (*AxInt)[0][j][k]=valX;
             (*AyInt)[0][j][k]=valY;
            SCDomainVal(Nx,j,k,0,1,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
            (*AyInt)[Nx][j][k]=valY;
            SCDomainVal(Nx-1,j,k,1,0,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
            (*AxInt)[Nx-1][j][k]=valX;
            }
       for (int k=0;k<=Nz-1;k++)//for (int k=0;k<=Nz;k+=sk)
          { int j=Ny;
           SCDomainVal(0,j,k,1,0,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
           (*AxInt)[0][j][k]=valX;
           (*AzInt)[0][j][k]=valZ;
           SCDomainVal(Nx,j,k,0,0,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
           (*AzInt)[Nx][j][k]=valZ;
           SCDomainVal(Nx-1,j,k,1,0,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
           (*AxInt)[Nx-1][j][k]=valX;
           }
            // plane paralel  XOZ (front and end planes)...........................................
            // we dont't repeat calculations on the ribs
         for (int i=0;i<=Nx-1;i+=si)
             for (int k=0;k<=Nz-1;k+=sk)
            {
            SCDomainVal(i,0,k,1,1,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
            (*AxInt)[i][0][k]=valX;
            (*AyInt)[i][0][k]=valY;
            (*AzInt)[i][0][k]=valZ;

            SCDomainVal(i,Ny,k,1,0,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
            (*AxInt)[i][Ny][k]=valX;
            (*AzInt)[i][Ny][k]=valZ;

            SCDomainVal(i,Ny-1,k,0,1,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
            (*AyInt)[i][Ny-1][k]=valY;
            }
//    additional line
         for (int i=0;i<=Nx-1;i++)
            {int k=Nz;
            SCDomainVal(i,0,k,1,1,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
            (*AxInt)[i][0][k]=valX;
            (*AyInt)[i][0][k]=valY;

            SCDomainVal(i,Ny,k,1,0,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
            (*AxInt)[i][Ny][k]=valX;

            SCDomainVal(i,Ny-1,k,0,1,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
            (*AyInt)[i][Ny-1][k]=valY;
            }

         for (int k=0;k<=Nz-1;k++)
        { int i=Nx;
        SCDomainVal(i,0,k,0,1,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
        (*AyInt)[i][0][k]=valY;
        (*AzInt)[i][0][k]=valZ;

        SCDomainVal(i,Ny,k,0,0,1,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
        (*AzInt)[i][Ny][k]=valZ;

        SCDomainVal(i,Ny-1,k,0,1,0,len,Jx,Jy,Jz,&valX,&valY,&valZ,n);
        (*AyInt)[i][Ny-1][k]=valY;
        }

//         char f_Ax[]="Ax_";
//         char f_Ay[]="Ay_";
//         char f_Az[]="Az_";

//         save_double_data3D(*AxInt,1003,0,f_Ax);
//         save_double_data3D(*AyInt,1003,0,f_Ay);
//         save_double_data3D(*AzInt,1003,0,f_Az);

// bilinear interpolation
//for each plane we foundother elements

//for Ax
 bilinearInterpXOY(AxInt,0,Nz,   Nx-1,Ny-1,si,sj);//
 bilinearInterpXOZ(AxInt,0,Ny,   Nx-1,Nz-1,si,sk);
 bilinearInterpYOZ(AxInt,0,Nx-1, Ny-1,Nz-1,sj,sk);


 //for Ay
  bilinearInterpXOY(AyInt,0,Nz,  Nx-1,Ny-1,si,sj);
  bilinearInterpXOZ(AyInt,0,Ny-1,Nx-1,Nz-1,si,sk);
  bilinearInterpYOZ(AyInt,0,Nx,  Ny-1,Nz-1,sj,sk);

  //for Az
   bilinearInterpXOY(AzInt,0,Nz-1, Nx-1,Ny-1,si,sj);
   bilinearInterpXOZ(AzInt,0,Ny,   Nx-1,Nz-1,si,sk);
   bilinearInterpYOZ(AzInt,0,Nx,   Ny-1,Nz-1,sj,sk);


//   save_double_data3D(*AxInt,1001,0,f_Ax);
//   save_double_data3D(*AyInt,1001,0,f_Ay);
//   save_double_data3D(*AzInt,1001,0,f_Az);

//// compare two solutions:
//  double ***AxIntE,***AyIntE,***AzIntE;

//  AxIntE=constrDouble3D();AyIntE=constrDouble3D();AzIntE=constrDouble3D();


//  IntegVectPotenInSpace(&AxIntE,&AyIntE,&AzIntE,len,Jx,Jy,Jz,n);

//  save_double_data3D(AxIntE,1000,0,f_Ax);
//  save_double_data3D(AyIntE,1000,0,f_Ay);
//  save_double_data3D(AzIntE,1000,0,f_Az);

//  double mdiff=0;
//  for (int i=0;i<Nx;i++)
//      for (int j=0;j<=Ny;j++)
//          if (mdiff<abs(AxIntE[i][j][0]-(*AxInt)[i][j][0]))
//          {
//              mdiff=abs(AxIntE[i][j][0]-(*AxInt)[i][j][0]);
//              cout<< i<<"  "<<j<<" is "<<AxIntE[i][j][0]<<"  and  "<<(*AxInt)[i][j][0]<<endl;
//          }

//  cout <<"mdif= "<<mdiff<<endl;
//  mdiff=0;
//  for (int i=0;i<Nx;i++)
//      for (int k=0;k<=Nz;k++)
//          if (mdiff<abs(AxIntE[i][0][k]-(*AxInt)[i][0][k]))
//          {
//              mdiff=abs(AxIntE[i][0][k]-(*AxInt)[i][0][k]);
//              cout<< i<<"  "<<k<<" is "<<AxIntE[i][0][k]<<"  and  "<<(*AxInt)[i][0][k]<<endl;
//          }

//  cout <<"mdif= "<<mdiff<<endl;
//  mdiff=0;
//  for (int k=0;k<=Nz;k++)
//      for (int j=0;j<=Ny;j++)
//          if (mdiff<abs(AxIntE[0][j][k]-(*AxInt)[0][j][k]))
//          {
//              mdiff=abs(AxIntE[0][j][k]-(*AxInt)[0][j][k]);
//              cout<< k<<"  "<<j<<" is "<<AxIntE[0][j][k]<<"  and  "<<(*AxInt)[0][j][k]<<endl;
//          }

//  cout <<"mdif= "<<mdiff<<endl;
////=====================Ay==================================
//  mdiff=0;
//  for (int i=0;i<=Nx;i++)
//      for (int j=0;j<Ny;j++)
//          if (mdiff<abs(AyIntE[i][j][0]-(*AyInt)[i][j][0]))
//          {
//              mdiff=abs(AyIntE[i][j][0]-(*AyInt)[i][j][0]);
//              cout<< i<<"  "<<j<<" is "<<AyIntE[i][j][0]<<"  and  "<<(*AyInt)[i][j][0]<<endl;
//          }

//  cout <<"mdif= "<<mdiff<<endl;
//  mdiff=0;
//  for (int i=0;i<=Nx;i++)
//      for (int k=0;k<Nz;k++)
//          if (mdiff<abs(AyIntE[i][0][k]-(*AyInt)[i][0][k]))
//          {
//              mdiff=abs(AyIntE[i][0][k]-(*AyInt)[i][0][k]);
//              cout<< i<<"  "<<k<<" is "<<AyIntE[i][0][k]<<"  and  "<<(*AyInt)[i][0][k]<<endl;
//          }

//  cout <<"mdif= "<<mdiff<<endl;
//  mdiff=0;
//  for (int k=0;k<=Nz;k++)
//      for (int j=0;j<Ny;j++)
//          if (mdiff<abs(AyIntE[0][j][k]-(*AyInt)[0][j][k]))
//          {
//              mdiff=abs(AyIntE[0][j][k]-(*AyInt)[0][j][k]);
//              cout<< k<<"  "<<j<<" is "<<AyIntE[0][j][k]<<"  and  "<<(*AyInt)[0][j][k]<<endl;
//          }

//  cout <<"mdif= "<<mdiff<<endl;
//  //================Az====================================
//  //=====================Ay==================================
//    mdiff=0;
//    for (int i=0;i<=Nx;i++)
//        for (int j=0;j<=Ny;j++)
//            if (mdiff<abs(AzIntE[i][j][0]-(*AzInt)[i][j][0]))
//            {
//                mdiff=abs(AzIntE[i][j][0]-(*AzInt)[i][j][0]);
//                cout<< i<<"  "<<j<<" is "<<AzIntE[i][j][0]<<"  and  "<<(*AzInt)[i][j][0]<<endl;
//            }

//    cout <<"mdif= "<<mdiff<<endl;
//    mdiff=0;
//    for (int i=0;i<=Nx;i++)
//        for (int k=0;k<Nz;k++)
//            if (mdiff<abs(AzIntE[i][0][k]-(*AzInt)[i][0][k]))
//            {
//                mdiff=abs(AzIntE[i][0][k]-(*AzInt)[i][0][k]);
//                cout<< i<<"  "<<k<<" is "<<AzIntE[i][0][k]<<"  and  "<<(*AzInt)[i][0][k]<<endl;
//            }

//    cout <<"mdif= "<<mdiff<<endl;
//    mdiff=0;
//    for (int k=0;k<Nz;k++)
//        for (int j=0;j<=Ny;j++)
//            if (mdiff<abs(AzIntE[0][j][k]-(*AzInt)[0][j][k]))
//            {
//                mdiff=abs(AzIntE[0][j][k]-(*AzInt)[0][j][k]);
//                cout<< k<<"  "<<j<<" is "<<AzIntE[0][j][k]<<"  and  "<<(*AzInt)[0][j][k]<<endl;
//            }

//    cout <<"mdif= "<<mdiff<<endl;


// freeDouble3D(AxIntE);freeDouble3D(AyIntE);freeDouble3D(AzIntE);

         double  end_time=clock();
         search_time=(double(end_time-start_time)/ (double)CLOCKS_PER_SEC);

         //cout<<"calculation of integral boundarz condition takes "<< search_time <<"   sec" <<endl;



}


 double Ax_ext_val(double y, double B0z, double a)
 {
     return (-0.5)*((B0z*y)*a);
 }

 double Ay_ext_val(double x,  double B0z, double a)
 {
     //double R = a / (2. * PI); //**** sleduet li zdes' zamenit' a na a+delta v sootvetstvii s formuloi na strochke 59
     //return ((-1.0)*B_ind*R*cos(s / R)); // cylinder
     return 0.5*B0z*(x*(2.-a)); // planar
 }



 void externalVectPoten(double ****AxExt,double ****AyExt,double ****AzExt, hh h, Cin c)
 {
     double xi,yj;
     double a=0.;// we can choose any value for a
     for (int i = 0; i < Nx + 1; i++)
                   for (int j = 0; j < Ny + 1; j++)
                       for (int k = 0; k < Nz + 1; k++)
                   {   xi= i*h.x;
                       yj= j*h.y;

                       (*AxExt)[i][j][k] = Ax_ext_val(yj, c.B0z,a);
                       (*AyExt)[i][j][k]=  Ay_ext_val(xi, c.B0z,a);
                       (*AzExt)[i][j][k] = 0.0;
                   }

     for (int i = 0; i < Nx + 1; i++)
                   for (int j = 0; j < Ny + 1; j++)
                       for (int k = 0; k < Nz + 1; k++)
                   {   xi= i*h.x;
                       yj= j*h.y;

                       (*AxExt)[i][j][k] += 0.;
                       (*AyExt)[i][j][k] +=  0.;
                       (*AzExt)[i][j][k] += c.B0x*yj;
                   }
     for (int i = 0; i < Nx + 1; i++)
                   for (int j = 0; j < Ny + 1; j++)
                       for (int k = 0; k < Nz + 1; k++)
                   {   xi= i*h.x;
                       yj= j*h.y;

                       (*AxExt)[i][j][k] += 0.;
                       (*AyExt)[i][j][k] += c.B0y*yj;
                       (*AzExt)[i][j][k] += 0.;
                   }
 }


void BT3D( double ****Bx,double ****By,double ****Bz,double ***Ax,double ***Ay,double ***Az,hh h,Cin c)
{// calculation of magnetic field
    //iDouble B[Nx+1][Ny+1];
    //B=ones(Nx,Ny)*c.B0+0i;
double ihx=1/h.x;
double ihy=1/h.y;
double ihz=1/h.z;


    fill3DReal(Bx,c.B0x);
    fill3DReal(By,c.B0y);
    fill3DReal(Bz,c.B0z);



    for (int i=0;i<=Nx-1;i++)
        for (int j=0;j<=Ny-1;j++)
            for (int k=0;k<=Nz-1;k++)
            {

                (*Bx)[i][j][k]=(Az[i][j+1][k]-Az[i][j][k])*ihy-(Ay[i][j][k+1]-Ay[i][j][k])*ihz;
                (*By)[i][j][k]=-(Az[i+1][j][k]-Az[i][j][k])*ihx+(Ax[i][j][k+1]-Ax[i][j][k])*ihz;
                (*Bz)[i][j][k]=(Ax[i][j][k]-Ax[i][j+1][k])*ihy+(Ay[i+1][j][k]-Ay[i][j][k])*ihx;
            }


 }

void DirichletCondX(double ****Ax,double ***AxInt )
{        //boundary condition  from the condition rot A=B0

    for (int j=0;j<=Ny;j++)
        for (int k=0;k<=Nz;k++)
            {
            (*Ax)[0][j][k]=AxInt[0][j][k];
            (*Ax)[Nx-1][j][k]=AxInt[Nx-1][j][k];
             }
    // plane paralel  XOZ
         for (int i=0;i<=Nx;i++)
             for (int k=0;k<=Nz;k++)
                {
                (*Ax)[i][0][k]=AxInt[i][0][k];
                (*Ax)[i][Ny][k]=AxInt[i][Ny][k];
                }

         for (int i=0;i<=Nx;i++)
             for (int j=0;j<=Ny;j++)
                {
                (*Ax)[i][j][0]=AxInt[i][j][0];
                (*Ax)[i][j][Nz]=AxInt[i][j][Nz];
                 }


}
void DirichletCondY(double ****Ay,double ***AyInt)
{        //boundary condition  from the condition rot A=B0

    for (int j=0;j<=Ny;j++)
        for (int k=0;k<=Nz;k++)
             {
            (*Ay)[0][j][k]=AyInt[0][j][k];
            (*Ay)[Nx][j][k]=AyInt[Nx][j][k];
             }

    // plane paralel  XOZ
         for (int i=0;i<=Nx;i++)
             for (int k=0;k<=Nz;k++)
                {
                (*Ay)[i][0][k]=AyInt[i][0][k];
                (*Ay)[i][Ny-1][k]=AyInt[i][Ny-1][k];
                }

         for (int i=0;i<=Nx;i++)
             for (int j=0;j<=Ny;j++)
                {
                (*Ay)[i][j][0]=AyInt[i][j][0];
                (*Ay)[i][j][Nz]=AyInt[i][j][Nz];
                 }

}

void DirichletCondZ(double ****Az,double ***AzInt )
{        //boundary condition  from the condition rot A=B0

    for (int j=0;j<=Ny;j++)
        for (int k=0;k<=Nz;k++)
            {
            (*Az)[0][j][k]=AzInt[0][j][k];
            (*Az)[Nx][j][k]=AzInt[Nx][j][k];
            }
    // plane paralel  XOZ
         for (int i=0;i<=Nx;i++)
             for (int k=0;k<=Nz;k++)
                {
                (*Az)[i][0][k]=AzInt[i][0][k];
                (*Az)[i][Ny][k]=AzInt[i][Ny][k];
                }

         for (int i=0;i<=Nx;i++)
             for (int j=0;j<=Ny;j++)
                {
                (*Az)[i][j][0]=AzInt[i][j][0];
                (*Az)[i][j][Nz-1]=AzInt[i][j][Nz-1];
                 }
}


void total_vect_potent(double ****Ax,double ****Ay,double ****Az,iDouble ****Ux,iDouble ****Uy,iDouble ****Uz,
                       double ***AxInd,double ***AyInd,double ***AzInd,double ***AxExt,double ***AyExt,double ***AzExt,Cin c, hh h)
{
    double khx=c.kappa*h.x;
    double khy=c.kappa*h.y;
    double khz=c.kappa*h.z;

for(int i=0;i<Nx+1;i++)
    for(int j=0;j<Ny+1;j++)
        for(int k=0;k<Nz+1;k++)
         {(*Ax)[i][j][k]=(AxInd)[i][j][k]+AxExt[i][j][k];
          (*Ay)[i][j][k]=(AyInd)[i][j][k]+AyExt[i][j][k];
          (*Az)[i][j][k]=(AzInd)[i][j][k]+AzExt[i][j][k];

          (*Ux)[i][j][k]=exp(-imI*(*Ax)[i][j][k]*khx);    //Ux=exp(-1i*Ax*c.kappa*h.x);
          (*Uy)[i][j][k]=exp(-imI*(*Ay)[i][j][k]*khy);    //Uy=exp(-1i*Ay*c.kappa*h.y);
          (*Uz)[i][j][k]=exp(-imI*(*Az)[i][j][k]*khz);
          }
}


void vector_potent(double ****AxInd,double ****AyInd,double ****AzInd,
                   double ***AxInt,double ***AyInt,double ***AzInt,
                   double ***Jx,double ***Jy,double ***Jz, hh h)

{   // we find only potential   indused current in superconductor

     double ***FAx,***FAy,***FAz;

     double inhx2,inhy2,inhz2;
     int flag_exit=0;

     double hmin2k=pow(min(h.x,min(h.y,h.z)),2.)/10.;

     inhx2=1/(h.x*h.x);
     inhy2=1/(h.y*h.y);
     inhz2=1/(h.z*h.z);

     FAx=constrDouble3D();
     FAy=constrDouble3D();
     FAz=constrDouble3D();

     int count_int=0;

     //   ........................Ax..........................................

    while (flag_exit==0)
// calculate operator FAx,FAy for input parameters  at each grid node............................................
        {
        for (int i=1;i<=Nx-2;i++)
            for (int j=1;j<=Ny-1;j++)
                for (int k=1;k<=Nz-1;k++)
                {
                   FAx[i][j][k]=Jx[i][j][k]+((*AxInd)[i+1][j][k]-2.0*(*AxInd)[i][j][k]+(*AxInd)[i-1][j][k])*inhx2
                                                +((*AxInd)[i][j+1][k]-2.0*(*AxInd)[i][j][k]+(*AxInd)[i][j-1][k])*inhy2
                                                +((*AxInd)[i][j][k+1]-2.0*(*AxInd)[i][j][k]+(*AxInd)[i][j][k-1])*inhz2;
                }
        sum3DReal((AxInd),FAx,h.tau);
        DirichletCondX(AxInd,  AxInt);


        //check exit condition
        if (max_abs_mass3D(FAx)<hmin2k)// && MaxAbs3D(*AxInd,Ax_old)<eps_B0)
            flag_exit=1;

        count_int++;
        }
    //cout<<"vector iter(x)="<<count_int<<"    discrepancy  "<<max_abs_mass3D(FAx)<<endl;

    //   ........................Ay..................................................................

    flag_exit=0;
    count_int=0;
    while (flag_exit==0)
        {
           for (int i=1;i<=Nx-1;i++)
               for (int j=1;j<=Ny-2;j++)
                   for (int k=1;k<=Nz-1;k++)
                   {
                        FAy[i][j][k]=Jy[i][j][k]+((*AyInd)[i+1][j][k]-2.0*(*AyInd)[i][j][k]+(*AyInd)[i-1][j][k])*inhx2
                                                +((*AyInd)[i][j+1][k]-2.0*(*AyInd)[i][j][k]+(*AyInd)[i][j-1][k])*inhy2
                                                +((*AyInd)[i][j][k+1]-2.0*(*AyInd)[i][j][k]+(*AyInd)[i][j][k-1])*inhz2;
                    }
            sum3DReal((AyInd),FAy,h.tau);
            //boundary and interface conditions
            DirichletCondY(AyInd,AyInt);

            //check exit condition
            if (max_abs_mass3D(FAy)<hmin2k)// && MaxAbs3D(*AyInd,Ay_old)<eps_B0)
                flag_exit=1;

            count_int++;
        }
    //cout<<"vector iter(z)="<<count_int<<"   discrepancy  "<<max_abs_mass3D(FAy)<<endl;

    //   ........................Az..................................................................
    flag_exit=0;
    count_int=0;

    while (flag_exit==0)
    {

          for (int i=1;i<=Nx-1;i++)
               for (int j=1;j<=Ny-1;j++)
                   for (int k=1;k<=Nz-2;k++)
                   {
                       FAz[i][j][k]=Jz[i][j][k]+((*AzInd)[i+1][j][k]-2.0*(*AzInd)[i][j][k]+(*AzInd)[i-1][j][k])*inhx2
                                                +((*AzInd)[i][j+1][k]-2.0*(*AzInd)[i][j][k]+(*AzInd)[i][j-1][k])*inhy2
                                                +((*AzInd)[i][j][k+1]-2.0*(*AzInd)[i][j][k]+(*AzInd)[i][j][k-1])*inhz2;

                   }
            sum3DReal((AzInd),FAz,h.tau);
            // boundary and interface conditions
            DirichletCondZ(AzInd,AzInt);

           // check exit condition(calculate relative error of residual)
            if (max_abs_mass3D(FAz)<hmin2k)
                flag_exit=1;

            count_int++;
    }
    //cout<<"vector iter(z)="<<count_int<<"   discrepancy  "<<max_abs_mass3D(FAz)<<endl;
//........................................................................................................
    freeDouble3D(FAx);
    freeDouble3D(FAy);
    freeDouble3D(FAz);
}
//====================================================================================================================================
// Получить количество тактов с момента последнего сброса процессора
static __inline__ unsigned long long rdtsc(void)
{
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

//====================================================================================================================================
void GL_VectorPoten3D()
{
    nn n; hh  h; Cin c;Coef acoef;
    koefGL kGL;

    int count=0,count_psi=0;
    int count_psi2=0;
    double tim=0;


    double ***Bx,***By,***Bz;
    double ***Ax,***Ay,***Az;
    double ***AxInd,***AyInd,***AzInd;
    double ***AxInt,***AyInt,***AzInt;
    double ***AxExt,***AyExt,***AzExt;

    double ***Jx, ***Jy,***Jz;
    double ***Jx_old, ***Jy_old,***Jz_old;

    iDouble ***Ux,***Uy,***Uz;

    iDouble ***psi;
    iDouble ***psi_old;                // psi_old=psi+ones(Nx+1,Ny+1);

    double ***phi;
    double ***phi_old;

    double ***len;

    double khx,khy,khz;


    int flag_psi=0;
    int flag_J=0;
    int del=0;
    double dis_psi=0;


    Bx=constrDouble3D();
    By=constrDouble3D();
    Bz=constrDouble3D();

    Ax=constrDouble3D();
    Ay=constrDouble3D();
    Az=constrDouble3D();

    AxExt=constrDouble3D();
    AyExt=constrDouble3D();
    AzExt=constrDouble3D();

    AxInt=constrDouble3D();
    AyInt=constrDouble3D();
    AzInt=constrDouble3D();

    AxInd=constrDouble3D();
    AyInd=constrDouble3D();
    AzInd=constrDouble3D();

    Ux=constrComp3D();
    Uy=constrComp3D();
    Uz=constrComp3D();

    psi=constrComp3D();

    Jx=constrDouble3D();
    Jy=constrDouble3D();
    Jz=constrDouble3D();
    Jx_old=constrDouble3D();
    Jy_old=constrDouble3D();
    Jz_old=constrDouble3D();

    psi=constrComp3D();
    psi_old=constrComp3D();       //psi_old=1; !!!!

    phi=constrDouble3D();
    phi_old=constrDouble3D();

    len=constrDouble3D();


    double *mas_err_psi;
    mas_err_psi=new double[Nt];


   //=======================================================================================


    //n=new nn[2];

    Initialise(&n,&h,&c,&acoef,&kGL);

    ostringstream dir;

    dir << "RC_Bz" << dlts(c.B0z)<<"_Bx"<< dlts(c.B0x)<<"_J"<<dlts(c.j_tr);

    int status;





    //status = mkdir(dir.str().c_str());//when we use Windows

    status = mkdir(dir.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);// when we use linux



    double hmin2=pow(min(h.x,min(h.y,h.z)),2.);
    double hmin=min(h.x,min(h.y,h.z));

    khx=c.kappa*h.x;
    khy=c.kappa*h.y;
    khz=c.kappa*h.z;


    CalculateLength(&len,h);

    save_general_inform3D(c,n,h,acoef,dir.str().c_str());

    externalVectPoten(&Ax,&Ay,&Az,h,c);


    srand(time(NULL));

        for (int i=n.sx;i<=n.ex;i++)
            for (int k=n.sz;k<=n.ez;k++)
                if(!((n.mx<i && i<=n.ex) && (n.mbz<k && k<n.muz)))
                    for (int j=n.sy;j<=n.ey;j++)
                             {
                                psi[i][j][k]=0.01*double(rand()%100)+imI*0.01*double(rand()%100);
                             }


        //count=10000;
        //char foldername[]="./R_Bz0_10_Bx0_00_J0_10/";

        //RecoverSolComp(psi,count,"psi_",foldername);
  //      RecoverSolReal(Ax,count,"Ax_",foldername);
  //      RecoverSolReal(Ay,count,"Ay_",foldername);
  //      RecoverSolReal(Az,count,"Az_",foldername);
  //      RecoverSolReal(phi,count,"phi_",foldername);
         // count++;flag_psi=1;


    for (int i = 0; i < Nx + 1; i++)
                  for (int j = 0; j < Ny + 1; j++)
                      for (int k = 0; k < Nz + 1; k++)
                  {
                      AxExt[i][j][k]=Ax[i][j][k];
                      AyExt[i][j][k]=Ay[i][j][k];
                      AzExt[i][j][k]=Az[i][j][k];

                      Ux[i][j][k]=exp(-imI*Ax[i][j][k]*khx);// we can calculate link variables only in sample !!!!
                      Uy[i][j][k]=exp(-imI*Ay[i][j][k]*khy);
                      Uz[i][j][k]=exp(-imI*Az[i][j][k]*khz);

                      psi_old[i][j][k]=psi[i][j][k]; //psi_old=psi;
                  }



    // ============================calculaion==============================




        while (count<Nt)
        {
//------------------------GL equation----------------------------------------------------
         GLeq3D( n,&psi,phi,Ux,Uy,Uz,h,c,&dis_psi);//        [ psi] = GLeq( n,psi,Ax,Ay,h,c);


         //if (count>9000 )//&& count%2==0)
         {
         mas_err_psi[count]=MaxAbs3D(psi,psi_old); //max(max(abs(psi-psi_old)));
         //cout<<"mas_err_psi["<<count<<"]="<< mas_err_psi[count]<<endl;

         if  (dis_psi<hmin2)//((mas_err_psi[count]<eps_psi))
            count_psi++;
         else
             count_psi=0;
         if (count_psi>=100)
                 flag_psi=1;//we found stable solution for vector potential
         // attune coefficient to recalculate GL time
         if (0)//(count>12500 && flag_psi==1)
             {
             double m_psi=max_mas_AbsComp3D(psi_old);
             if  ((mas_err_psi[count]/m_psi<0.01))
                 count_psi2++;
             else
                 count_psi2=0;

              if (count_psi2>100)
              {
                  kGL.koef=min(kGL.koef+1,kGL.koef_max);
                  cout<<"coefficient of Ginzburg Landau time ="<<kGL.koef<<endl;
                  h.t=kGL.koef*hmin2;
                  count_psi2=0;
              }
               if  ((mas_err_psi[count]/m_psi>0.2))
               {
                  kGL.koef=max(kGL.koef-1,kGL.koef_min);
                  h.t=kGL.koef*hmin2;
               }
             }
            }



             for(int i=n.sx;i<=n.ex;i++)
                 for(int k=n.sz;k<=n.ez;k++)
                         if(!((n.mx<i && i<=n.ex) && (n.mbz<k && k<n.muz)))
                             for(int j=n.sy;j<=n.ey;j++)
                                psi_old[i][j][k]=psi[i][j][k]; //psi_old=psi;

//            if (count%10000==0)
//                    cout<<"iteration = "<<count<<endl;

//----------------------Vect Potenc -------------------------------------------------------
             if ((count>Niter) && (flag_psi==1)) //   condition should provide stable solution of GL(number more than Niter and value at the last 100 iterations differs none more than eps.B)
                //if   (count % 1==0)
                {
                     sc_current(&Jx,&Jy,&Jz, psi,Ux,Uy,Uz,c, n, h);


                     //flag_J=MaxAbs3DReal3Arr(Jx,Jx_old,Jy,Jy_old,Jz,Jz_old,hmin2);

                     if (flag_J==0 || count%50==0)
                     {
                         IntegVectPotenInSpaceMod( &AxInt,&AyInt,&AzInt,len,Jx,Jy,Jz,n);

                         for(int i=0;i<Nx+1;i++)
                             for(int j=0;j<Ny+1;j++)
                                 for(int k=0;k<Nz+1;k++)
                                         {
                                         Jx_old[i][j][k]=Jx[i][j][k];
                                         Jy_old[i][j][k]=Jy[i][j][k];
                                         Jz_old[i][j][k]=Jz[i][j][k];
                                         }
                     }
                     //================================================================================================
                        //double time1 = clock() / (double)CLOCKS_PER_SEC;

                      vector_potent(&AxInd,&AyInd,&AzInd,AxInt,AyInt,AzInt,Jx,Jy,Jz,h); //find Ax,Ay

                      //double time2 = clock() / (double)CLOCKS_PER_SEC;
                      //double cpu_time = time2 - time1;      // Потраченное процессорное время
                      //<<"CPU TIME: "<<cpu_time<<endl;
                      //cout<<"vector potential, count="<<count<<endl;
                    //}

 // -----------------Scalar potencial-------------------------------------------------

                    //if ((count>Niter) && (flag_psi==1)) //  we include scalar potential when we find stable solution of GL+vect potential (iteration number more than Niter and value at the last 100 iterations differs none more than eps.B
                    //if   (count % 10==0)
                    //  {
                     if(del==0)
                     {
                     char f_psi[]="psi_";
                     save_complex_data3D(psi,count,tim,f_psi,dir.str().c_str());
                     del=1;
                     }
                     scalar_poten3D(n,psi,&phi,Ux,Uy,Uz, h,c,1);

                     //cout<<"scalar potential, count="<<count<<endl;
                      //  }

                     // total vector potential and linc variables
                    total_vect_potent(&Ax,&Ay,&Az,&Ux,&Uy,&Uz,AxInd,AyInd,AzInd,AxExt,AyExt,AzExt,c,h);
         }


        if ((count>=10000) && (count%100==0)&& (flag_psi==1))
            {// save basic calculated matrix
            BT3D( &Bx,&By,&Bz,Ax,Ay,Az,h,c);

            char f_abs_psi[]="abs_psi_";
            char f_Bx[]="Bx_";
            char f_By[]="By_";
            char f_Bz[]="Bz_";

            save_abs_data3D(psi,count,tim,f_abs_psi,dir.str().c_str());
            save_double_data3D(Bx,count,tim,f_Bx,dir.str().c_str());
            save_double_data3D(By,count,tim,f_By,dir.str().c_str());
            save_double_data3D(Bz,count,tim,f_Bz,dir.str().c_str());

            char f_phi[]="phi_";


            save_double_data3D(phi,count,tim,f_phi,dir.str().c_str());
            char fvolt[]="voltage";
            double sum1=0.,sum2=0.;
            for (int j=n.sy;j<=n.ey;j++)
                {
                    for (int k=n.sz;k<=n.mbz;k++)
                    {
                   sum1+=phi[n.ex][j][k];
                    }
                    for (int k=n.muz;k<=n.ez;k++)
                    {
                    sum2+=phi[n.ex][j][k];
                    }
                }
            sum1*=1./((double)(n.ey-n.sy)*(n.mbz-n.sz));
            sum2*=1./((double)(n.ey-n.sy)*(n.ez-n.muz));

            save_voltage3D( (sum1-sum2),count,fvolt,dir.str().c_str());

            char fmas_dis[]="mas_discr";

            save_voltage3D( dis_psi,count,fmas_dis,dir.str().c_str());

            char f_Ax[]="Ax_";
            char f_Ay[]="Ay_";
            char f_Az[]="Az_";

            save_double_data3D(Ax,count,tim,f_Ax,dir.str().c_str());
            save_double_data3D(Ay,count,tim,f_Ay,dir.str().c_str());
            save_double_data3D(Az,count,tim,f_Az,dir.str().c_str());

            char f_psi[]="psi_";
            char f_Jx[]="Jx_";
            char f_Jy[]="Jy_";
            char f_Jz[]="Jz_";

            save_complex_data3D(psi,count,tim,f_psi,dir.str().c_str());

            save_double_data3D(Jx,count,tim,f_Jx,dir.str().c_str());
            save_double_data3D(Jy,count,tim,f_Jy,dir.str().c_str());
            save_double_data3D(Jz,count,tim,f_Jz,dir.str().c_str());


            }
         tim=tim+h.t;
         count++;
            }

//---------------------------- free -------------------------------------------------------------------
        delete[] mas_err_psi;


        freeDouble3D(Bx);
        freeDouble3D(By);
        freeDouble3D(Bz);

        freeDouble3D(Ax);
        freeDouble3D(Ay);
        freeDouble3D(Az);

        freeDouble3D(AxExt);
        freeDouble3D(AyExt);
        freeDouble3D(AzExt);


        freeDouble3D(AxInt);
        freeDouble3D(AyInt);
        freeDouble3D(AzInt);

        freeDouble3D(AxInd);
        freeDouble3D(AyInd);
        freeDouble3D(AzInd);

        freeComp3D(Ux);
        freeComp3D(Uy);
        freeComp3D(Uz);


        freeDouble3D(Jx);
        freeDouble3D(Jy);
        freeDouble3D(Jz);

        freeComp3D(psi);
        freeComp3D(psi_old);

        freeDouble3D(phi);
        freeDouble3D(phi_old);
        freeDouble3D(len);
}
//====================================================================================================================================
//====================================================================================================================================
