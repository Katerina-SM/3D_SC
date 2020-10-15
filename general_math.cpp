#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <time.h>
#include <cstring>

#include <string>

#include <algorithm>

#include <sys/types.h>

#include <sys/stat.h>



#include "general_math.h"

//using namespace std;


// Получить количество тактов с момента последнего сброса процессора
static __inline__ unsigned long long rdtsc(void)
{
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}



void sum(iDouble ***A, iDouble **B,double k)
{


    for(int i=0;i<Nx+1;i++)
        for(int j=0;j<Ny+1;j++)
            (*A)[i][j]=(*A)[i][j]+k*B[i][j];


}
void sum3D(iDouble ****A, iDouble ***B,double koef)
{


    for(int i=0;i<Nx+1;i++)
        for(int j=0;j<Ny+1;j++)
            for(int k=0;k<Nz+1;k++)
            {

                (*A)[i][j][k]=(*A)[i][j][k]+koef*B[i][j][k];

            }



}

void sum3DReal(double ****A, double ***B,double koef)
{


    for(int i=0;i<Nx+1;i++)
        for(int j=0;j<Ny+1;j++)
            for(int k=0;k<Nz+1;k++)
            {

                (*A)[i][j][k]=(*A)[i][j][k]+koef*B[i][j][k];

            }



}
int isNull3D(iDouble ***B)
{


    for(int i=0;i<Nx+1;i++)
        for(int j=0;j<Ny+1;j++)
            for(int k=0;k<Nz+1;k++)
                if (B[i][j][k]!=0.)
                           return 0;// array B contain npnzero elements!!!! something wrong!!

    return 1;// all elements of array B is 0



}

void fill(iDouble ***B,double val)
{


    for(int i=0;i<Nx+1;i++)
        for(int j=0;j<Ny+1;j++)
            (*B)[i][j]=val;


}

void fill3D(iDouble ****B,double val)
{


    for(int i=0;i<Nx+1;i++)
        for(int j=0;j<Ny+1;j++)
            for(int k=0;k<Nz+1;k++)
                (*B)[i][j][k]=val;


}
void fill3DReal(double ****B,double val)
{


    for(int i=0;i<Nx+1;i++)
        for(int j=0;j<Ny+1;j++)
            for(int k=0;k<Nz+1;k++)
                (*B)[i][j][k]=val;


}

void abs_val(iDouble **B,double ***abs_B)
{


    for(int i=0;i<Nx+1;i++)
        for(int j=0;j<Ny+1;j++)
            (*abs_B)[i][j]=abs(B[i][j]);


}
void fun_tt(char tt[8],int count)
{
    //tt=new char[8];

//definition of process number in char-type; for the file-name
int i_1000000, i_100000, i_10000, i_1000, i_100, i_10, i_1;
i_1000000 = count/1000000;
i_100000 = (count-1000000*i_1000000)/100000;
i_10000 = (count-1000000*i_1000000-100000*i_100000)/10000;
i_1000 = (count-1000000*i_1000000-100000*i_100000-10000*i_10000)/1000;
i_100 = (count-1000000*i_1000000-100000*i_100000-10000*i_10000-1000*i_1000)/100;
i_10  = (count-1000000*i_1000000-100000*i_100000-10000*i_10000-1000*i_1000-100*i_100)/10;
i_1   = (count-1000000*i_1000000-100000*i_100000-10000*i_10000-1000*i_1000-100*i_100-10*i_10);
tt[0] = i_1000000 + 48;
tt[1] = i_100000  + 48;
tt[2] = i_10000   + 48;
tt[3] = i_1000    + 48;
tt[4] = i_100     + 48;
tt[5] = i_10      + 48;
tt[6] = i_1       + 48;
tt[7] = '\0';
}


double max_mass(double **mas, int count_s, int count_y)
{
    double max = mas[0][0];

    for (int i = 0; i < count_s; i++)
    {
        for (int j = 0; j < count_y; j++)
        {
            if (mas[i][j] > max)
            {
                max = mas[i][j];
            }
        }
    }
    return max;
}
double max_mass3D(double ***mas)
{
    double max = mas[0][0][0];

    for (int i = 0; i <= Nx; i++)
    {
        for (int j = 0; j <= Ny; j++)

            for (int k = 0; k <= Nz; k++)
            {
                if (mas[i][j][k] > max)
                {
                    max = mas[i][j][k];
                }
            }
    }
    return max;
}
double max_abs_mass3D(double ***mas)
{
    double max = mas[0][0][0];

    for (int i = 0; i <= Nx; i++)
    {
        for (int j = 0; j <= Ny; j++)

            for (int k = 0; k <= Nz; k++)
            {   double buf=abs(mas[i][j][k]);
                if ( buf  > max)
                {
                    max = buf;
                }
            }
    }
    return max;
}

double max_mas_AbsComp(iDouble **mas)
{
    double max = abs(mas[0][0]);

    for (int i = 0; i <= Nx; i++)
    {
        for (int j = 0; j <= Ny; j++)
        {
            if (abs(mas[i][j]) > max)
            {
                max = abs(mas[i][j]);
            }
        }
    }
    return max;
}
double max_mas_AbsComp3D(iDouble ***mas)
{
    double max = abs(mas[0][0][0]);

    for (int i = 0; i <= Nx; i++)
    {
        for (int j = 0; j <= Ny; j++)

            for (int k = 0; k <= Nz; k++)
            {
                if (abs(mas[i][j][k]) > max)
                {
                    max = abs(mas[i][j][k]);
                }
            }
    }
    return max;
}

void save_voltage3D(double var,int count,char fname[],string dir)
{
    char *out_psi;
    out_psi=new char[50];

    //ostringstream directoryCurr;
    //directoryCurr<<"Results";



    ostringstream fullPath;


    strcpy(out_psi,fname);
    strcat(out_psi,".txt\0");



    fullPath.str("");
    fullPath<<dir<<"/"<<out_psi;

    //ofstream out00(out_psi,ios::app);
    ofstream out00(fullPath.str().c_str(),ios::app);
    
    if(!out00)
    cerr<<"Error at the opening!!!.\n"<<endl;


    out00 << var << '\t'<<count<<'\t' << endl;

    out00.close();
   delete []out_psi;
}



 void save_abs_data(iDouble **var,int count,double time,char fname[],string dir)
 {
     char *out_psi;
     out_psi=new char[50];

     char tt[8];
     //definition of process number in char-type; for the file-name
      fun_tt(tt,count);

     //ostringstream directoryCurr;
     //directoryCurr<<"Results";

     //int status=mkdir(directoryCurr.str().c_str());

     ostringstream fullPath;
//     //fullPath<<directoryCurr.str()<<"/"<<coeff_out;

     //ofstream out(fullPath.str().c_str(),ios::app);



     strcpy(out_psi,fname);//"mod_PSI_\0");
     //strcpy(out_psi,itoa(count,buff,10));
     strcat(out_psi,tt);//itoa(count,buff,10));
     strcat(out_psi,".txt\0");

     fullPath.str("");
     //fullPath<<directoryCurr.str()<<"/"<<out_psi;
     fullPath<<dir<<"/"<<out_psi;

     ofstream out;//(fullPath.str().c_str(),ios::app);
     out.open(fullPath.str().c_str());

    if(!out)
        cerr<<"Error at the opening!!!.\n"<<endl;
    out << time*tau*pow(10,9) << "  " << endl;

    for (int i = 0; i < Nx + 1; i++)
    {
        for (int j = 0; j < Ny + 1; j++)
        {
        out <<  abs(var[i][j]) << '\t';
            if(j != Ny)
                out << "  ";
            else
                out << endl;
        }
    }
    out << endl <<endl;
    out.close();
    delete []out_psi;
 }

 void save_abs_data3D(iDouble ***var,int count,double time,char fname[],string dir)
 {
     char *out_psi;
     out_psi=new char[50];

     char tt[8];
     //definition of process number in char-type; for the file-name
      fun_tt(tt,count);

     //ostringstream directoryCurr;
     //directoryCurr<<"Results";

     //int status=mkdir(directoryCurr.str().c_str());

     ostringstream fullPath;
//     //fullPath<<directoryCurr.str()<<"/"<<coeff_out;

     //ofstream out(fullPath.str().c_str(),ios::app);



     strcpy(out_psi,fname);//"mod_PSI_\0");
     //strcpy(out_psi,itoa(count,buff,10));
     strcat(out_psi,tt);//itoa(count,buff,10));
     strcat(out_psi,".txt\0");

     fullPath.str("");
     //fullPath<<directoryCurr.str()<<"/"<<out_psi;
     fullPath<<dir<<"/"<<out_psi;

     ofstream out;//(fullPath.str().c_str(),ios::app);
     out.open(fullPath.str().c_str());

    if(!out)
        cerr<<"Error at the opening!!!.\n"<<endl;
    out << time*tau*pow(10,9) << "  " << endl;

// save k matrix of Nx*Ny one by one
for (int k = 0; k < Nz + 1; k++)
{
    for (int i = 0; i < Nx + 1; i++)
    {
        for (int j = 0; j < Ny + 1; j++)

        {
        out <<  abs(var[i][j][k]) << '\t';
            if(j != Ny )//|| k!=Nz)
                out << "  ";
            else
                out << endl;
        }
    }
}

    out << endl <<endl;
    out.close();
    delete []out_psi;
 }

 void save_double_data(double **var,int count,double time,char fname[],string dir)
 {
     char *out_psi;
     out_psi=new char[50];

     char tt[8];
     //definition of process number in char-type; for the file-name
     fun_tt(tt,count);

     //ostringstream directoryCurr;
     //directoryCurr<<"Results";

     //int status=mkdir(directoryCurr.str().c_str());

     ostringstream fullPath;
//     //fullPath<<directoryCurr.str()<<"/"<<coeff_out;

     //ofstream out00(fullPath.str().c_str(),ios::app);


     strcpy(out_psi,fname);
     strcat(out_psi,tt);
     strcat(out_psi,".txt\0");



     fullPath.str("");
     fullPath<<dir<<"/"<<out_psi;

     //ofstream out00(out_psi,ios::app);
     ofstream out00;//(fullPath.str().c_str(),ios::app);
     out00.open(fullPath.str().c_str());

     if(!out00)
     cerr<<"Error at the opening!!!.\n"<<endl;


     out00 << time*tau*pow(10,9) << "  " << endl;

     for (int i = 0; i < Nx + 1; i++)
     {
         for (int j = 0; j < Ny + 1; j++)
         {
         out00 <<  var[i][j] << '\t';
             if(j != Ny)
                 out00 << "  ";
             else
                 out00 << endl;
         }
     }

     out00 << endl << endl;
     out00.close();
    delete []out_psi;
 }

 void save_double_data3D(double ***var,int count,double time,char fname[],string dir)
 {
     char *out_psi;
     out_psi=new char[50];

     char tt[8];
     //definition of process number in char-type; for the file-name
     fun_tt(tt,count);

     //ostringstream directoryCurr;
     //directoryCurr<<"Results";

     //int status=mkdir(directoryCurr.str().c_str());

     ostringstream fullPath;
//     //fullPath<<directoryCurr.str()<<"/"<<coeff_out;

     //ofstream out00(fullPath.str().c_str(),ios::app);


     strcpy(out_psi,fname);
     strcat(out_psi,tt);
     strcat(out_psi,".txt\0");



     fullPath.str("");
     fullPath<<dir<<"/"<<out_psi;

     //ofstream out00(out_psi,ios::app);
     ofstream out00;//(fullPath.str().c_str(),ios::app);
     out00.open(fullPath.str().c_str());

     if(!out00)
     cerr<<"Error at the opening!!!.\n"<<endl;


     out00 << time*tau*pow(10,9) << "  " << endl;
     // save k matrix of Nx*Ny one by one
  for (int k = 0; k < Nz + 1; k++)
     for (int i = 0; i < Nx + 1; i++)
     {
         for (int j = 0; j < Ny + 1; j++)
         {
         out00 <<  var[i][j][k] << '\t';
             if(j != Ny)
                 out00 << "  ";
             else
                 out00 << endl;
         }
     }

     out00 << endl << endl;
     out00.close();
    delete []out_psi;
 }


 void save_general_inform3D(Cin c,nn n,hh h,Coef acoef, string  dir)
 {
     char *out_psi;
     out_psi=new char[50];
     char fname[]="general_inf";

     ostringstream fullPath;


     strcpy(out_psi,fname);
     //strcat(out_psi,tt);
     strcat(out_psi,".txt\0");



     fullPath.str("");

     fullPath<<dir<<"/"<<out_psi;

     ofstream out00;
     out00.open(fullPath.str().c_str());
     if(!out00)
     cerr<<"Error at the opening!!!.\n"<<endl;


     out00  << "nondimensional domain        \t"<<ax<<" * "<<ay<<" * "<<az << '\t'<<endl;
     out00  << "dimensional domain size(nm)  \t"<<ax*lambda*pow(10,9)<<" * "<<ay*lambda*pow(10,9)<<" * "<<az*lambda*pow(10,9) <<endl;
     out00  << " ============================================================================== "<<endl;
     out00  << "grid nodes is                  \t"<<Nx<<" * "<<Ny<<" * "<<Nz<<'\t'<<endl;
     out00  << "step size is  h.x="<<h.x<<", h.y="<<h.y<<", h.z="<<h.z<<", h.t="<<h.t<<", h.tau="<<h.tau<<'\t'<<endl;
     out00  << "dimensional step size (nm)   h.x="<<h.x*lambda*pow(10,9)<<", h.y="<<h.y*lambda*pow(10,9)<<", h.z="<<h.z*lambda*pow(10,9)<<", h.t="<<h.t*lambda*pow(10,9)<<", h.tau="<<h.tau*lambda*pow(10,9)<<'\t'<<endl;
     out00  << " ============================================================================== "<<endl;
     out00  << "nondimensional size of sample convex hull   \t"<<(n.ex-n.sx)*ax/Nx<<" * "<<(n.ey-n.sy)*ay/Ny<<" * "<<(n.ez-n.sz)*az/Nz<< '\t'<<endl;
     out00  << "dimensional size of sample convex hull(nm)  \t"<<(n.ex-n.sx)*ax/Nx*lambda*pow(10,9)<<" * "<<(n.ey-n.sy)*ay/Ny*lambda*pow(10,9)<<" * "<<(n.ez-n.sz)*az/Nz*lambda*pow(10,9)<< '\t'<<endl;
     out00  <<"n.sx= "<<n.sx << ", n.mx= "<<n.mx<<",n.ex= "<<n.ex<<endl;
     out00  <<"n.sy= "<<n.sy <<", n.ey= "<<n.ey<<endl;
     out00  <<"n.sz= "<<n.sz << ", n.mbz= "<<n.mbz<<",n.mub= "<<n.muz<<",n.ez= "<<n.ez<<endl;

     out00 << " ============================================================================== "<<endl;

     out00  << "nondimensional magnetic field  B0x \t"<<c.B0x<< '\t'<<endl;
     out00  << "dimensional magnetic field(mT) B0x \t"<<c.B0x*sqrt(2)*Hc<< '\t'<<endl;

     out00  << "nondimensional magnetic field  B0y \t"<<c.B0y<< '\t'<<endl;
     out00  << "dimensional magnetic field(mT) B0y \t"<<c.B0y*sqrt(2)*Hc<< '\t'<<endl;

     out00  << "nondimensional magnetic field  B0z \t"<<c.B0z<< '\t'<<endl;
     out00  << "dimensional magnetic field(mT) B0z \t"<<c.B0z*sqrt(2)*Hc<< '\t'<<endl;
     out00  << " ============================================================================== "<<endl;

     out00  << "coherence length (nm) \t"<<xi*pow(10,9)<< '\t'<<endl;

     out00  << "GL parameter kappa    \t"<<c.kappa<< '\t'<<endl;
     out00  << "penetration depth (nm)\t"<<lambda*pow(10,9)<< '\t'<<endl;
     out00  << "vortex size (nm)      \t"<<2*xi*pow(10,9)<< '\t'<<endl;

      out00  << " ============================================================================== "<<endl;

     out00  << "nondimensional conductivity coefficient            \t"<<c.sigma<< '\t'<<endl;
     out00  << "dimensional conductivity coefficient (1/(muOm*m))  \t"<<c.sigma*sigma_0*pow(10,-6)<< '\t'<<endl;
      out00  << " ============================================================================== "<<endl;

     out00  << "nondimensional current density       \t"<<c.j_tr<< '\t'<<endl;
     out00  << "dimensional current density(GA/m^2)  \t"<<c.j_tr*j_0*pow(10,-9)<< '\t'<<endl;
     out00  << " ============================================================================== "<<endl;
     out00  << " axternal vector potential has the form (-1/2*Bz(y-sy)*a;1/2*Bz(x-sx)*(2-a);0) "<<endl;
     out00  << " a="<<acoef.a<<" , sx= "<<acoef.sx<<", sy="<<acoef.sy<<endl;


     out00 << endl << endl;
     out00.close();
    delete []out_psi;
 }





 void save_complex_data(iDouble **var,int count,double time,char fname[],string dir)
 {
     char *out_psi;
     out_psi=new char[50];

     //double abs_psi[Nx+1][Ny+1];

     //char buff[15];
     char tt[8];
     //definition of process number in char-type; for the file-name
      fun_tt(tt,count);

     //ostringstream directoryCurr;
     //directoryCurr<<"Results";

     //int status=mkdir(directoryCurr.str().c_str());

     ostringstream fullPath;
//     //fullPath<<directoryCurr.str()<<"/"<<coeff_out;

     //ofstream out(fullPath.str().c_str(),ios::app);



     strcpy(out_psi,fname);//"PSI_\0");
     strcat(out_psi, tt);
     strcat(out_psi,".txt\0");

     fullPath.str("");
     //fullPath<<directoryCurr.str()<<"/"<<out_psi;
     fullPath<<dir<<"/"<<out_psi;

     ofstream out;//(fullPath.str().c_str(),ios::app);
     out.open(fullPath.str().c_str());
    if(!out)
        cerr<<"Error at the opening!!!.\n"<<endl;

    out << time*tau*pow(10,9) << "  " << endl;

    for (int i = 0; i < Nx + 1; i++)
    {
        for (int j = 0; j < Ny + 1; j++)
        {
        out <<  real(var[i][j])<<"  "<<imag(var[i][j]) << '\t';
            if(j != Ny)
                out << "  ";
            else
                out << endl;
        }
    }



    out << endl <<endl;
    out.close();
    delete []out_psi;
 }

 void save_complex_data3D(iDouble ***var,int count,double time,char fname[],string dir)
 {
     char *out_psi;
     out_psi=new char[50];

     //double abs_psi[Nx+1][Ny+1];

     //char buff[15];
     char tt[8];
     //definition of process number in char-type; for the file-name
      fun_tt(tt,count);

//     ostringstream directoryCurr;
//     directoryCurr<<"Results";

     //int status=mkdir(directoryCurr.str().c_str());

     ostringstream fullPath;
//     //fullPath<<directoryCurr.str()<<"/"<<coeff_out;

     //ofstream out(fullPath.str().c_str(),ios::app);



     strcpy(out_psi,fname);//"PSI_\0");
     strcat(out_psi, tt);
     strcat(out_psi,".txt\0");

     fullPath.str("");
//     fullPath<<directoryCurr.str()<<"/"<<out_psi;
     fullPath<<dir<<"/"<<out_psi;

     ofstream out;//(fullPath.str().c_str(),ios::app);
     out.open(fullPath.str().c_str());

    if(!out)
        cerr<<"Error at the opening!!!.\n"<<endl;

    out << time*tau*pow(10,9) << "  " << endl;

 // save k matrix of Nx*Ny one by one
 for (int k = 0; k < Nz + 1; k++)
    for (int i = 0; i < Nx + 1; i++)
    {
        for (int j = 0; j < Ny + 1; j++)
        {
        out <<  real(var[i][j][k])<<"  "<<imag(var[i][j][k]) << '\t';
            if(j != Ny)
                out << "  ";
            else
                out << endl;
        }
    }



    out << endl <<endl;
    out.close();
    delete []out_psi;
 }

 iDouble ***read_complex3D(char fname[])
 {
     double a;
     double ra,ia;


     iDouble ***var;
     var=constrComp3D();


      ifstream inread;
      if(!inread)
        {cout<<"Error of opening file!!!\n";
        return var;}

      inread.open(fname,ios::in );
      if (inread.is_open())
      {


        inread >>a;
        for (int k=0;k<=Nz;k++)
          for (int i=0;i<=Nx;i++)
              for (int j=0;j<=Ny;j++)
              {
                  inread >>ra;
                  inread >>ia;
                  (var)[i][j][k]=ra+imI*ia;

              }

        inread.close();
      }

      else cout << "Unable to open file";
      return var;

 }

 double ***read_double3D(char fname[])
 {
     double a;



     double ***var;
     var=constrDouble3D();


      ifstream inread;
      if(!inread)
        {cout<<"Error of opening file!!!\n";
        return var;}

      inread.open(fname,ios::in );
      if (inread.is_open())
      {


        inread >>a;
        for (int k=0;k<=Nz;k++)
          for (int i=0;i<=Nx;i++)
              for (int j=0;j<=Ny;j++)
              {

                  inread >>(var)[i][j][k];

              }

        inread.close();
      }

      else cout << "Unable to open file";
      return var;

 }

 iDouble **read_complex(char fname[])
 {
     double a;
     double ra,ia;


     iDouble **var= new iDouble*[Nx + 1];
     for (int i = 0; i < Nx + 1; i++)
         (var)[i] = new iDouble[Ny + 1];


      ifstream inread;
      if(!inread)
        {cout<<"Error of opening file!!!\n";
        return var;}

      inread.open(fname,ios::in );
      if (inread.is_open())
      {


          inread >>a;
          for (int i=0;i<=Nx;i++)
              for (int j=0;j<=Ny;j++)
              {
                  inread >>ra;
                  inread >>ia;
                  (var)[i][j]=ra+imI*ia;

              }

        inread.close();
      }

      else cout << "Unable to open file";
      return var;

 }


 double **read_double(char fname[])
 {
     double a;



     double **var= new double*[Nx + 1];
     for (int i = 0; i < Nx + 1; i++)
         (var)[i] = new double[Ny + 1];


      ifstream inread;
      if(!inread)
        {cout<<"Error of opening file!!!\n";
        return var;}

      inread.open(fname,ios::in );
      if (inread.is_open())
      {


          inread >>a;
          for (int i=0;i<=Nx;i++)
              for (int j=0;j<=Ny;j++)
              {
                  inread >>(var)[i][j];

              }

        inread.close();
      }

      else cout << "Unable to open file";
      return var;

 }

// void save_test()
// {
//     char *out_psi;
//     out_psi=new char[50];
//     double p;


//     ostringstream directoryCurr;
//     directoryCurr<<"Results";

//     ostringstream fullPath;

//     strcpy(out_psi,"example_.txt\0");

//     fullPath.str("");
//     fullPath<<directoryCurr.str()<<"/"<<out_psi;

//     ofstream out(fullPath.str());
//    if(!out)
//        cerr<<"Error at the opening (double) !!!.\n"<<endl;
//    p=0.34567;
//    out << p << "  ";

//    for (int i=0;i<=5;i++)
//        for (int j=0;j<=5;j++)
//        {
//            p=double((i+j)/10.0);
//            out <<p<<'\t';
//            if (j!=5)
//                    out<<"  ";
//            else
//            out <<endl;
//        }
//    out << endl <<endl;
//    out.close();
//    delete []out_psi;
// }

 double MaxAbs(iDouble **A,iDouble **B)
 {
     double res=abs(A[0][0]-B[0][0]);

     for(int i=0;i<Nx+1;i++)
         for(int j=0;j<Ny+1;j++)
             if (abs(A[i][j]-B[i][j])>res)
                 res=abs(A[i][j]-B[i][j]);

     return res;
 }
 double MaxAbs3D(iDouble ***A,iDouble ***B)
 {
     double res=abs(A[0][0][0]-B[0][0][0]);

     for(int i=0;i<Nx+1;i++)
         for(int j=0;j<Ny+1;j++)
             for(int k=0;k<Nz+1;k++)
             if (abs(A[i][j][k]-B[i][j][k])>res)
             {
                 res=abs(A[i][j][k]-B[i][j][k]);

             }

     return res;
 }
 double MaxAbs3DReal(double ***A,double ***B)
 {
     double res=abs(A[0][0][0]-B[0][0][0]);

     for(int i=0;i<Nx+1;i++)
         for(int j=0;j<Ny+1;j++)
             for(int k=0;k<Nz+1;k++)
             if (abs(A[i][j][k]-B[i][j][k])>res)
             {
                 res=abs(A[i][j][k]-B[i][j][k]);

             }

     return res;
 }

 int MaxAbs3DReal3Arr(double ***A,double ***A_old,double ***B,double ***B_old,double ***C,double ***C_old, double hmin2)
 {
     int flag_J=0;
     double resA=abs(A[0][0][0]-A_old[0][0][0]);
     double resB=abs(B[0][0][0]-B_old[0][0][0]);
     double resC=abs(C[0][0][0]-C_old[0][0][0]);

     for(int i=0;i<Nx+1;i++)
         for(int j=0;j<Ny+1;j++)
             for(int k=0;k<Nz+1;k++)
             {
             if (abs(A[i][j][k]-A_old[i][j][k])>resA)
             {
                 resA=abs(A[i][j][k]-A_old[i][j][k]);
              }
             if (abs(B[i][j][k]-B_old[i][j][k])>resB)
             {
                 resB=abs(B[i][j][k]-B_old[i][j][k]);
              }
             if (abs(C[i][j][k]-C_old[i][j][k])>resC)
             {
                 resC=abs(C[i][j][k]-C_old[i][j][k]);
              }

            }
     flag_J=(resA<hmin2 && resB<hmin2 && resC<hmin2);
     return flag_J;
 }

void copyComp3D(iDouble ***toCopy,iDouble ****Copy)
{

    for(int i=0;i<Nx+1;i++)
        for(int j=0;j<Ny+1;j++)
            for(int k=0;k<Nz+1;k++)
            {
                (*Copy)[i][j][k]=toCopy[i][j][k];
            }


}



 iDouble **constrComp()
 {

    iDouble **var= new iDouble*[Nx + 1];
    for (int i = 0; i < Nx + 1; i++)
        (var)[i] = new iDouble[Ny + 1];

    for (int i=0;i<Nx+1;i++)
        for (int j=0;j<Ny+1;j++)
             var[i][j]=0.0;

    return var;

 }
 void freeComp(iDouble **var)
 {

      for (int i = 0; i < Nx + 1; i++)
          delete[] var[i];
      delete[] var;

 }
 void freeDouble(double **var)
{

     for (int i = 0; i < Nx + 1; i++)
         delete[] var[i];
     delete[] var;
 }

 double **constrDouble()
 {

    double **var= new double*[Nx + 1];
    for (int i = 0; i < Nx + 1; i++)
        (var)[i] = new double[Ny + 1];

    for (int i=0;i<Nx+1;i++)
        for (int j=0;j<Ny+1;j++)
             var[i][j]=0.0;

    return var;
 }

 double ***constrDouble3D()
 {

//    double ***var;
//    var=new double**[Nx + 1];

//    for (int i = 0; i < Nx + 1; i++)
//    {
//        (var)[i] = new double*[Ny + 1];
//        for (int j = 0; j < Ny + 1; j++)
//            (var)[i][j] = new double[Nz + 1];
//    }


//    for (int i=0;i<Nx+1;i++)
//        for (int j=0;j<Ny+1;j++)
//            for (int k=0;k<Nz+1;k++)
//             var[i][j][k]=0.0;


    double ***var= (double ***)calloc(Nx+1, sizeof(double **));

    for (int i = 0; i < Nx + 1; i++)
    {
        (var)[i] = (double**)calloc(Ny+1,sizeof(double **));
        for (int j = 0; j < Ny + 1; j++)
            (var)[i][j] = (double *)calloc(Nz+1, sizeof(double));
    }




    return var;
 }


 void freeDouble3D(double ***var)
{

//     for (int i = 0; i < Nx + 1; i++)
//     {
//         for (int j = 0; j < Ny + 1; j++)
//            delete[] var[i][j];
//         delete[] var[i];
//     }

//     delete[] var;

          for (int i = 0; i < Nx + 1; i++)
          {
              for (int j = 0; j < Ny + 1; j++)
                 free(var[i][j]);
              free(var[i]);
          }

          free(var);
 }

 iDouble ***constrComp3D()
 {

//    iDouble ***var;
//    var=new iDouble**[Nx + 1];

//    for (int i = 0; i < Nx + 1; i++)
//    {
//        (var)[i] = new iDouble*[Ny + 1];
//        for (int j = 0; j < Ny + 1; j++)
//            (var)[i][j] = new iDouble[Nz + 1];
//    }


//    for (int i=0;i<Nx+1;i++)
//        for (int j=0;j<Ny+1;j++)
//            for (int k=0;k<Nz+1;k++)
//             var[i][j][k]=0.0;


// we allocate memory and set var[i][j][k]=0! using  calloc function
    iDouble ***var= (iDouble ***)calloc(Nx+1, sizeof(iDouble **));

    for (int i = 0; i < Nx + 1; i++)
    {
        (var)[i] = (iDouble**)calloc(Ny+1,sizeof(iDouble **));
        for (int j = 0; j < Ny + 1; j++)
            (var)[i][j] = (iDouble *)calloc((Nz+1), sizeof(iDouble));
    }

    return var;
 }



 void freeComp3D(iDouble ***var)
{

//     for (int i = 0; i < Nx + 1; i++)
//     {
//         for (int j = 0; j < Ny + 1; j++)
//            delete[] var[i][j];
//         delete[] var[i];
//     }

//     delete[] var;
     for (int i = 0; i < Nx + 1; i++)
     {
         for (int j = 0; j < Ny + 1; j++)
            free(var[i][j]);
         free(var[i]);
     }

     free(var);
 }


 iDouble As_value(double s, double y, double z, double B_ind)
 {
     return (0.0);
 }

 iDouble Ay_value(double s, double y, double z, double B_ind, double a)
 {
     //double R = a / (2. * PI); //**** sleduet li zdes' zamenit' a na a+delta v sootvetstvii s formuloi na strochke 59
     //return ((-1.0)*B_ind*R*cos(s / R)); // cylinder
     return ((-1.0)*B_ind*s); // planar
 }

 iDouble Az_value(double s, double y, double z, double B_ind)
 {
     return (0.0);
 }



 double ReDss(double psiRE_nest, double psiIM_nest, double psi_mid, double psiRE_prev, double psiIM_prev, double ReUs_mid, double ImUs_mid, double ReUs_prev, double ImUs_prev, double hs)
 {
     return ((psiRE_nest*ReUs_mid + psiIM_nest*ImUs_mid - 2.0*psi_mid + psiRE_prev*ReUs_prev - psiIM_prev*ImUs_prev) / (hs*hs));
 }

 double ImDss(double psiRE_nest, double psiIM_nest, double psi_mid, double psiRE_prev, double psiIM_prev, double ReUs_mid, double ImUs_mid, double ReUs_prev, double ImUs_prev, double hs)
 {
     return (((-1.0)*psiRE_nest*ImUs_mid + psiIM_nest*ReUs_mid - 2.0*psi_mid + psiIM_prev*ReUs_prev + psiRE_prev*ImUs_prev) / (hs*hs));
 }

 double ReDyy(double psiRE_nest, double psiIM_nest, double psi_mid, double psiRE_prev, double psiIM_prev, double ReUy_mid, double ImUy_mid, double ReUy_prev, double ImUy_prev, double hy)
 {
     return ((psiRE_nest*ReUy_mid + psiIM_nest*ImUy_mid - 2.0*psi_mid + psiRE_prev*ReUy_prev - psiIM_prev*ImUy_prev) / (hy*hy));
 }

 double ImDyy(double psiRE_nest, double psiIM_nest, double psi_mid, double psiRE_prev, double psiIM_prev, double ReUy_mid, double ImUy_mid, double ReUy_prev, double ImUy_prev, double hy)
 {
     return (((-1.0)*psiRE_nest*ImUy_mid + psiIM_nest*ReUy_mid - 2.0*psi_mid + psiIM_prev*ReUy_prev + psiRE_prev*ImUy_prev) / (hy*hy));
 }
 string dlts (double x)
 {

 string s;
 s=to_string(int(x))+'_';
 if ((int)(x*10)==0)
     s+='0';
 s+=to_string((int)((x-int(x))*100));
 return s;
}

 void RecoverSolComp(iDouble ***X,int count,char fname[],char out_psi[] )
 {
 char tt[8];


 fun_tt(tt,count);

 strcat(out_psi,fname);
 strcat(out_psi,tt);
 strcat(out_psi,".txt\0");


 X=read_complex3D(out_psi);

 }
 void RecoverSolReal(double ***X,int count,char fname[],char out_psi[] )
 {
 char tt[8];


 fun_tt(tt,count);

 strcat(out_psi,fname);
 strcat(out_psi,tt);
 strcat(out_psi,".txt\0");


 X=read_double3D(out_psi);

 }
