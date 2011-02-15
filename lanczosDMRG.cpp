/**
 * @file lanczosDMRG.cpp
 *
 * @brief Lanczos routine for performing ED in DMRG
 *
 * Roger Melko Aug. 4 2006
 */
#include <cmath>
#include <iomanip>
#include "blitz/array.h"
#include "exceptions.h"
#include "lanczosDMRG_helpers.h"
#include "tqli2.h"
#include "matrixManipulation.h"
#include "lanczosDMRG.h"

BZ_USING_NAMESPACE(blitz)

/**
 * @brief A function to reduce the Hamiltonian to a tri-diagonal form  
 * 
 */
int LanczosED(blitz::Array<double,2>& Ham, blitz::Array<double,1>& Psi, double *En, const int N)
{
  int MAXiter, EViter;
  int min;
  int Lexit;
  double E0;

  int STARTIT=3; //iteration which diagonz. begins
  int LIT=100;   //max number of Lanczos iterations
  
  //Matrices
  blitz::Array<double,1> V0(N);  
  blitz::Array<double,1> Vorig(N);
  blitz::Array<double,1> V1(N);  //Ground state vector
  blitz::Array<double,1> V2(N);
  blitz::Array<double,1> alpha(LIT);
  blitz::Array<double,1> beta(LIT);
  //For ED of tri-di Matrix routine (C)
  int nn, rtn;
  blitz::Array<double,1> e(LIT);
  blitz::Array<double,1> d(LIT); 
  blitz::Array<double,2> Hmatrix(LIT,LIT);

  //tensor indices
  firstIndex i;    secondIndex j;
  thirdIndex k;   
  
  int iter = 0;
  //
  // initialize with randon numbers are normalize
  //
  randomize(Vorig);
  Normalize(Vorig);  

  for (EViter = 0; EViter < 2; EViter++) {//0=get E0 converge, 1=get eigenvec

    iter = 0;
    V0 = Vorig;
    
    if (EViter == 1) Psi = V0*(Hmatrix(0,min));
    
    V1 = 0;
    beta(0)=0;  //beta_0 not defined
    
    V1 = sum(Ham(i,j)*V0(j),j); // V1 = H |V0> 
    
    alpha(0) = dotProduct(V0,V1);
    
    V1 -= alpha(0)*V0;
    beta(1) = calculateNorm(V1);

    if (fabs(pow(beta(1),2)) < 0.000000001){   //wavefnt Ham alread GS
      *En = alpha(0);
      return 1;
    }

    V1 /= beta(1);

    if (EViter == 1) Psi += V1*(Hmatrix(1,min));
    
    // done 0th iteration
    
    Lexit = 0;   //exit flag
    E0 = 1.0;    //previous iteration GS eigenvalue
    while(Lexit != 1){
      
      iter++;
      
      V2 = sum(Ham(i,j)*V1(j),j); // V2 = H |V1>
      //V2 -= beta(iter)*V0;
      
      alpha(iter) = dotProduct(V1,V2);
      
      V2 = V2-alpha(iter)*V1 -  beta(iter)*V0;
      beta(iter+1) = calculateNorm(V2);
      
      V2 /=beta(iter+1);

      if (EViter == 1) {Psi += V2*(Hmatrix(iter+1,min));
	//cout<<Psi<<" S \n";
      }
      
      V0 = V1;
      V1 = V2;
      
      if (iter > STARTIT && EViter == 0){
	
	//diagonalize tri-di matrix
	d(0) = alpha(0);
	for (int ii=1;ii<=iter;ii++){
	  d(ii) = alpha(ii);
	  e(ii-1) = beta(ii);
	}
	e(iter) = 0;
	
	nn = iter+1;
	rtn = tqli2(d,e,nn,Hmatrix,0);
	
	min = 0;
	for (int ii=1;ii<=iter;ii++)
	  if (d(ii) < d(min))  min = ii;
	
	if ( (E0 - d(min)) < 1E-5) {
	  Lexit = 1;
// 	  cout<<"Lanc :"<<iter<<" ";
// 	  cout<<setprecision(12)<<d(min)<<"\n";  
	  *En = d(min);
	}
	else {
	  E0 = d(min);
	}

        if (iter == LIT-2) {
          LIT += 100;
          cout<<LIT<<" Resize Lan. it \n";
          d.resize(LIT);
          e.resize(LIT);
          Hmatrix.resize(LIT,LIT);
          alpha.resizeAndPreserve(LIT);
          beta.resizeAndPreserve(LIT);
        }//end resize
	
      }//end STARTIT

      if (EViter == 1 && iter == MAXiter) Lexit = 1;
      
    }//while
   
    if (EViter == 0){
      MAXiter = iter;
      //diagonalize tri-di matrix
      d(0) = alpha(0);
      for (int ii=1;ii<=iter;ii++){
	d(ii) = alpha(ii);
	e(ii-1) = beta(ii);
      }
      e(iter) = 0;
      //calculate eigenvector
      Hmatrix = 0;
      for (int ii=0;ii<=iter;ii++)
	Hmatrix(ii,ii) = 1.0; //identity matrix
      nn = iter+1;
      rtn = tqli2(d,e,nn,Hmatrix,1);
      min = 0;
      for (int ii=1;ii<=iter;ii++)
	if (d(ii) < d(min))  min = ii;
    }
    
  }//repeat (EViter) to transfrom eigenvalues H basis
  
  Normalize(Psi);

//   cout<<Psi<<" Psi \n";

//   V2 = sum(Ham(i,j)*Psi(j),j);
//   for (int ii=0;ii<N;ii++)
//     cout<<ii<<" "<<V2(ii)/Psi(ii)<<" EVdiv \n";

  return 0;
} 
/**
 * @brief A function to calculate the ground state function using the
 * Lanczos algorithm
 *
 * @param Hm is a 4-index tensor with the Hamiltonian
 * @param Ed is a matrix with the result of the calculation
 * @param nn is ??
 *
 * Returns the ground state eigenvalue and eigenvector using the 
 * Lanczos function
 */
double calculateGroundState(Array<double,4>& Hm, Array<double,2>& Ed)
{
    const int nn=sqrt(Hm.numElements());

    Array<double,2> Ham2d(nn,nn);

    //complicated integer square root?
    int L = static_cast<int>(std::sqrt(1.0*nn));          

    Ham2d=reduceM2M2(Hm);

    Array<double,1> Psi(nn);  //return eigenvector
    double En;                //return eigenvalue

    int lrt = LanczosED(Ham2d, Psi, &En, nn); 
    if (lrt == 1) 
      throw dmrg::Exception("Lanczos early term error");

    //repack Psi as 2D Matrix - Eigenvector
    int c2 = 0;
    for (int i1=0; i1<L; i1++)
    {
      for (int i2=0; i2<L; i2++)
      {
	  double melem = Psi(c2);
	  //if (melem < 1e-16 && melem > -1e-16) melem = 0;
	  Ed(i2,i1) = melem;
	  c2++;
      }
    }
    return En;  //ground state eigenvalue
}
//end lanczosDMRG.cpp
