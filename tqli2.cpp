/**
 * @file tqli2.cpp
 *
 * @brief Implementation for the tqli2 function using blitz arrays
 *
 * @author Roger Melko 
 * @author Ivan Gonzalez
 * @date $Date$
 *
 * $Revision$ 
 */
#include <cmath>
#include <iomanip>
#include "blitz/array.h"
#include "tqli2.h"

#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
/**
 * @brief A function to diagonalize a tridiagonal matrix
 *
 * @param d an array with the elements in the diagonal of the matrix
 * @param e an array with the elements in the off-diagonal of the matrix
 * @param n an int with the size of the array d
 * @param z an matrix with the eigenvectors of the diagonal matrix
 * @param Evects an int to control whether you calculate the eigenvectors
 * or not
 *
 * @return an int with code to check success (1) or failure (0). (Note the
 * evil use of the flag.)
 *
 * Diagonalizes a tridiagonal matrix: d[] is input as the diagonal elements,
 * e[] as the off-diagonal.  If the eigenvalues of the tridiagonal matrix
 * are wanted, input z as the identity matrix.  If the eigenvalues of the
 * original matrix reduced by tred2 are desired, input z as the matrix
 * output by tred2.  The kth column of z returns the normalized eigenvectors,
 * corresponding to the eigenvalues output in d[k].
 * April 2005, Roger Melko, modified from Numerical Recipies in C v.2
 * Modified from www.df.unipi.it/~moruzzi/
 * Feb 23 2005: modified to use Blitz++ arrays
 */
int tqli2(blitz::Array<double,1>& d, blitz::Array<double,1>& e, int n, 
	blitz::Array<double,2>& z, const int Evects)
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
	    std::cout <<"Too many iterations in tqli() \n";
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
//end tqli2.cpp
