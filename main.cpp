#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <time.h>
#include <ctime>
#include <cstring>

#include <complex>
//#include <tgmath.h>
#include <vector>
#include <algorithm>

#include "gl_vectorpoten.h"



#include <sys/types.h>
#include <sys/stat.h>

#include <omp.h>




//void fun(double A[5][5])
//{
//    for(int i=0;i<5;i++)
//        for(int j=0;j<5;j++)
//            A[i][j]=i+j;

//}

//void fun2(double ***BB)
//{
//    (*BB)=(double**)malloc(sizeof(double*)*5);
//    for(int i=0;i<5;i++)
//        (*BB)[i]=(double*)malloc(sizeof(double)*5);


//    for(int i=0;i<5;i++)
//        for(int j=0;j<5;j++)
//            (*BB)[i][j]=i+j;

//}

//void fun22(double **BB)
//{
//    (*BB)=(double*)malloc(sizeof(double)*5);
//    for(int i=0;i<5;i++)
//        //for(int j=0;j<5;j++)
//            (*BB)[i]=i+5;

//}
//void test()
//{
//    double A[9][5]={10};
//    double *BB;  // dinamic array
//    double **AA; // dinamic matrix


//    fun(A);// return static array

//    fun2(&AA);//return dinamic matrix

//    fun22(&BB);//return dinamic array

////  iDouble x(0, 1);
////  iDouble x2(1, -2);
////  iDouble z;

////  z=abs(x2);


////  //z=x*x2*imI;
////  std::cout<<z<<std::endl;
////  iDouble  a[Nx+1][Ny+1];
////  fill(a,0.1);

//  //std::cout<<A<<std::endl;



//}



//double  exampleOMP()
//{
//    double a=0, b=1000;
//    int n=10000000;

//     double h = (b - a) / n;
//     double sum = 0.0;
//     #pragma omp parallel
//     {
//         double sumloc = 0.0;
//         #pragma omp for
//         for (int i = 0; i < n; i++)
//             sumloc += pow(a + h * (i + 0.5),2.);

//         #pragma omp atomic
//         sum += sumloc;
//     }
//     sum *= h;
//     return sum;

// }
//double  example()
//{
//    double a=0, b=1000;
//    int n=10000000;

//    double h = (b - a) / n;
//    double sum = 0.0;
//    for (int i = 0; i < n; i++)
//        sum += pow(a + h * (i + 0.5),2.);
//    sum *= h;
//    return sum;
//}

int main(int argc, char *argv[])
{
//    double t1=0,t2=0;
//    for (int i=1;i<100;i++)
//        {
//        double t = omp_get_wtime();
//        exampleOMP();
//        t = omp_get_wtime() - t;        //printf("Thread %d execution time: %.6f sec.\n",               omp_get_thread_num(), t);
//        cout<< "time is "<< t<<endl;

//        t1+=t;

//        t = omp_get_wtime();
//        example();
//        t = omp_get_wtime() - t;        //printf("Thread %d execution time: %.6f sec.\n",               omp_get_thread_num(), t);
//        cout<< "time without OMP is "<< t<<endl;
//        t2+=t;
//        }
//cout<<"average time is "<<t1/100<<endl;
//cout<<"average time without OMP is "<<t2/100<<endl;




 //=================================================================================================================

GL_VectorPoten3D();

 //=================================================================================================================


 //system("pause");

 return 0;
}


