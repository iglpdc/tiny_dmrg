//Oct. 26, Roger Melko ORNL
//
//Exact diagonalization program for a real, symmetric matrix
//Householder reduction to tridiagonal form, then ED
//
//Takes input as Blitz++ Arrays<>
//
//set the below value to "1" if you want the eigenvectors, and to "0" if not
#define EVECTS 1  

#include <blitz/array.h>
#include "tred3.h"
#include "tqli2.h"
#include "HHtridi8.h"

/**
 * @brief A function calculate the density matrix eigenvalues
 *
 * Returns the M largest eigenvectors 
 *
 */
void DMlargeEigen(blitz::Array<double,2>& Hm, blitz::Array<double,2>& Od, const int nn, 
		  const int mm)
{
  blitz::Array<double,1> e(nn);
  blitz::Array<double,1> d(nn); 

  int L = nn/2;

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
  
  int i,j;
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
   for (int kk=0; kk<mm; kk++){
     Addx += d(j);
     j--;
   }
//   std::cout<<"Addx "<<" "<<1.0-Addx<<endl;

  j=nn-1;
  for (int  kk=0; kk<mm; kk++){
//      std::cout<<d[j]<<" | ";
    for (i=0; i<nn; i++){
      Od(kk,i) = Hm(i,inx(j));   //define the truncation matrix
//        std::cout<<Hmatrix[i][inx[j]]<<" ";
    }
    j--;
//      std::cout<<endl;
  } 

}//EigenValues of density matrix
