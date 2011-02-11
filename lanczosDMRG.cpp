/**
 * @file lanczosDMRG.cpp
 *
 * @brief Lanczos routine for performing ED in DMRG
 *
 * Roger Melko Aug. 4 2006
 */
#include <cmath>
#include <iomanip>
#include "lanczosDMRG.h"

/**
 * @brief A function to calculate the ground state function using the
 * Lanczos algorithm
 *
 * @param Hm is a 4-index tensor with the Hamiltonian
 * @param Ed is a matrix with the result of the calculation
 * @param nn is ??
 *
 * Returns the ground state eigenvalue and eigenvector using the Lanczos function
 *
 */
double calculateGroundState(Array<double,4>& Hm, Array<double,2>& Ed, const int nn)
{
  Array<double,2> Ham2d(nn,nn);

  //complicated integer square root?
  int L = static_cast<int>(std::sqrt(1.0*nn));          

  int c1=0;
  int c2;

  for (int i1=0; i1<L; i1++)
  {
    for (int i2=0; i2<L; i2++)
    {
      c2=0;
      for (int i3=0; i3<L; i3++)
      {
	for (int i4=0; i4<L; i4++)
	{
	  Ham2d(c1,c2) = Hm(i1,i2,i3,i4);  //pack as 2D matrix
	  c2++;
	}
      }
      c1++;
    }
  }
  
  Array<double,1> Psi(nn);  //return eigenvector
  double En;                //return eigenvalue

  int lrt = LanczosED(Ham2d, Psi, &En, nn); 
  if (lrt == 1) cout<<" Lanczos early term error \n)";

  //repack Psi as 2D Matrix - Eigenvector
  c2 = 0;
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

/**
 * @brief A function to reduce the Hamiltonian to a tri-diagonal form  
 * 
 */
int LanczosED(Array<double,2>& Ham, Array<double,1>& Psi, double *En, const int N)
{
  int MAXiter, EViter;
  int min;
  int Lexit;
  double E0;

  int STARTIT=3; //iteration which diagonz. begins
  int LIT=100;   //max number of Lanczos iterations
  
  //Matrices
  Array<double,1> V0(N);  
  Array<double,1> Vorig(N);
  Array<double,1> V1(N);  //Ground state vector
  Array<double,1> V2(N);
  Array<double,1> alpha(LIT);
  Array<double,1> beta(LIT);
  //For ED of tri-di Matrix routine (C)
  int nn, rtn;
  Array<double,1> e(LIT);
  Array<double,1> d(LIT); 
  Array<double,2> Hmatrix(LIT,LIT);

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


#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
/**
 * @brief A function to diagonalize a tridiagonal matrix
 *
 * April 2005, Roger Melko, modified from Numerical Recipies in C v.2
 * modified from www.df.unipi.it/~moruzzi/
 * Diagonalizes a tridiagonal matrix: d[] is input as the diagonal elements,
 * e[] as the off-diagonal.  If the eigenvalues of the tridiagonal matrix
 * are wanted, input z as the identity matrix.  If the eigenvalues of the
 * original matrix reduced by tred2 are desired, input z as the matrix
 * output by tred2.  The kth column of z returns the normalized eigenvectors,
 * corresponding to the eigenvalues output in d[k].
 * Feb 23 2005: modified to use Blitz++ arrays
 */
int tqli2(Array<double,1>& d, Array<double,1>& e, int n, Array<double,2>& z, 
	const int Evects)
{
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;

  for (l=0;l<n;l++) {
    iter=0;
    do { 
      for (m=l;m<n-1;m++) { 
	dd=fabs(d(m))+fabs(d(m+1));
	if (fabs(e(m))+dd == dd) break;
      }
      if (m!=l) { 
	if (iter++ == 30) { 
	  cout <<"Too many iterations in tqli() \n";
	  return 0;
	}
	g=(d(l+1)-d(l))/(2.0*e(l));
	r=sqrt((g*g)+1.0);
	g=d(m)-d(l)+e(l)/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) { 
	  f=s*e(i);
	  b=c*e(i);
	  if (fabs(f) >= fabs(g)) { 
	    c=g/f;r=sqrt((c*c)+1.0);
	    e(i+1)=f*r;
	    c *= (s=1.0/r);
	  }
	  else { 
	    s=f/g;r=sqrt((s*s)+1.0);
	    e(i+1)=g*r;
	    s *= (c=1.0/r);
	  }
	  g=d(i+1)-p;
	  r=(d(i)-g)*s+2.0*c*b;
	  p=s*r;
	  d(i+1)=g+p;
	  g=c*r-b;
	  /*EVECTS*/
	  if (Evects == 1) {
	    for (k=0;k<n;k++) { 
	      f=z(k,i+1);
	      z(k,i+1)=s*z(k,i)+c*f;
	      z(k,i)=c*z(k,i)-s*f;
	    }
	  }//Evects
	}
	d(l)=d(l)-p;
	e(l)=g;
	e(m)=0.0;
      }
    } while (m!=l);
  }
  return 1;
}
//end lanczosDMRG.cpp
