//Oct. 26, Roger Melko ORNL
//
//Exact diagonalization program for a real, symmetric matrix
//Householder reduction to tridiagonal form, then ED
//
//Takes input as Blitz++ Arrays<>
//
//set the below value to "1" if you want the eigenvectors, and to "0" if not
#define EVECTS 1  

#include <cmath>
#include <blitz/array.h>
#include "tred3.h"
#include "tqli2.h"
#include "HHtridi8.h"

//BZ_USING_NAMESPACE(blitz)
//using namespace std;

/******************************************************************/
void DMlargeEigen(blitz::Array<double,2>& Hm, blitz::Array<double,2>& Od, const int nn, 
		  const int mm)
/*******************************************************************
 *** returns the M largest eigenvectors ***************************/
{
  int i,j,kk;
  int L;
    
  blitz::Array<double,1> e(nn);
  blitz::Array<double,1> d(nn); 

  L = nn/2;

  //Householder reduction
  tred3(Hm,d,e,nn);

  //Iterative diagonalization
  int rtn = tqli2(d,e,nn,Hm,1);

  //now Hmatrix[j][i] contains the eigenvectors corresponding
  //to d[i]

  // output eigenvalues d[] and corresponding evectors
  // m LARGEST

   double Esum = 0;
   for (int i1=0; i1<nn; i1++){
     Esum += d(i1);
//      cout<<d[i1]<<" | ";
//      for (int i2=0; i2<nn; i2++)
//        cout<<Hmatrix[i2][i1]<<" ";
//      cout<<endl;
   }
   if (Esum < 0.999999 || Esum > 1.00001)
     std::cout<<"Esum error:"<<Esum<<std::endl;	

//  int *inx = new int[nn];
  blitz::Array<int,1> inx(nn);
  for (j=0; j<nn; j++) inx(j)=j;
  double a;
  int b;
  for (j=1; j<nn; j++){         //STRIAGHT INSERTION SORT O(N^2)    
    a = d(j);
    b = inx(j);
    i=j-1;
    while (i>=0 && d(i) >a){
      d(i+1)=d(i);
      inx(i+1)=inx(i);
      i--;
    }
    d(i+1)=a;
    inx(i+1)=b;
  }

   double Addx = 0;
   j=nn-1;
   for (kk=0; kk<mm; kk++){
     Addx += d(j);
     j--;
   }
//   std::cout<<"Addx "<<" "<<1.0-Addx<<endl;

  j=nn-1;
  for (kk=0; kk<mm; kk++){
//      std::cout<<d[j]<<" | ";
    for (i=0; i<nn; i++){
      Od(kk,i) = Hm(i,inx(j));   //define the truncation matrix
//        std::cout<<Hmatrix[i][inx[j]]<<" ";
    }
    j--;
//      std::cout<<endl;
  } 

}//EigenValues of density matrix

///*********************************************************************/ 
/////
///// @brief Householder reduces a real symmetric matrix a to tridiagonal form
/////
///// @param a is the matrix
/////
///// On output, a is replaced by the orthogonal matrix effecting the transformation.
///// The diagonal elements are stored in d[], and the offdiagonal elements are
///// stored in e[].
///// March 30 2006: modified to use Blitz++ arrays
/////
//void tred3(Array<double,2>& a, Array<double,1>& d, Array<double,1>& e, int n)
//{
//  int i,j,k,l;
//  double f,g,h,hh,scale;
// 
//  for (i=n-1;i>0;i--)
//    { l=i-1;h=scale=0.0;
//      if (l>0)
//        { for (k=0;k<=l;k++) scale+=fabs(a(i,k));
//          if (scale==0.0) {e(i)=a(i,l);continue;}// .....skip transformation
//          // .............................. use scaled a's for transformation
//          for (k=0;k<=l;k++) { a(i,k)/=scale;h+=a(i,k)*a(i,k);}
//          f=a(i,l);
//          g=sqrt(h);
//          if (f>0.0) g=-g;
//          e(i)=scale*g;
//          h-=f*g;
//          a(i,l)=f-g;
//          f=0.0;
//          for (j=0;j<=l;j++)
//            { a(j,i)=a(i,j)/h;g=0.0;
//              for (k=0;k<=j;k++) g+=a(j,k)*a(i,k);
//              for (k=j+1;k<=l;k++) g+=a(k,j)*a(i,k);
//              e(j)=g/h;
//              f += e(j)*a(i,j);
//            }
//          hh=f/(h+h);
//          for (j=0;j<=l;j++)
//            { f=a(i,j);e(j)=g=e(j)-hh*f;
//              for (k=0;k<=j;k++) a(j,k)-=(f*e(k)+g*a(i,k));
//            }
//        }
//      else e(i)=a(i,l);
//      d(i)=h;
//    }
//                          // eigenvectors
//#if (EVECTS == 1)
//  d(0)=0.0;e(0)=0.0;
//  for (i=0;i<n;i++)
//    { l=i-1;
//      if (d(i))
//        { for (j=0;j<=l;j++)
//            { g=0.0;
//              for (k=0;k<=l;k++) g+=a(i,k)*a(k,j);
//              for (k=0;k<=l;k++) a(k,j)-=g*a(k,i);
//            }
//        }
//      d(i)=a(i,i);
//      a(i,i)=1.0;
//      for (j=0;j<=l;j++) a(j,i)=a(i,j)=0.0;
//    }
//#else 
//  for (i=0;i<n;i++) d(i)=a(i,i);
//#endif
//
//  for (i=0;i<n-1;i++) e(i)=e(i+1);
//  e(n-1)=0;
//}


