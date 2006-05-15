//May 2005, R. Melko UCSB
//
//Exact diagonalization program for a real, symmetric matrix
//Householder reduction to tridiagonal form, then ED

//set the below value to "1" if you want the eigenvectors, and to "0" if not
#define EVECTS 1  

using namespace std;
#include <iostream>
#include <fstream>
#include <cmath>

void EigenValues(double **Hmatrix, double *d, const int nn)
{
  int i,j;
  
  void tred2(double **, double *, double *,int );
  int tqli(double *, double *, int, double **);
  
  double *e = new double[nn];  

  //Householder reduction
  tred2(Hmatrix,d,e,nn);

  //Iterative diagonalization
  int rtn = tqli(d,e,nn,Hmatrix);

  //now Hmatrix[j][i] contains the eigenvectors corresponding
  //to d[i]

  delete [] e;  

}//EigenValues

/*********************************************************************/ 
 
void tred2(double **a,double *d,double *e,int n)
/***
 May 2005, Roger Melko, modified from Numerical Recipies in C v.2
 modified from www.df.unipi.it/~moruzzi/
 Householder reduces a real symmetric matrix **a to tridiagonal form.
 On output, a is replaced by the orthogonal matrix effecting the transformation. 
 The diagonal elements are stored in d[], and the offdiagonal elements are
 stored in e[].
***/
{
  int i,j,k,l;
  double f,g,h,hh,scale;
 
  for (i=n-1;i>0;i--)
    { l=i-1;h=scale=0.0;
      if (l>0)
        { for (k=0;k<=l;k++) scale+=fabs(a[i][k]);
          if (scale==0.0) {e[i]=a[i][l];continue;}// .....skip transformation
          // .............................. use scaled a's for transformation
          for (k=0;k<=l;k++) { a[i][k]/=scale;h+=a[i][k]*a[i][k];}
          f=a[i][l];
          g=sqrt(h);
          if (f>0.0) g=-g;
          e[i]=scale*g;
          h-=f*g;
          a[i][l]=f-g;
          f=0.0;
          for (j=0;j<=l;j++)
            { a[j][i]=a[i][j]/h;g=0.0;
              for (k=0;k<=j;k++) g+=a[j][k]*a[i][k];
              for (k=j+1;k<=l;k++) g+=a[k][j]*a[i][k];
              e[j]=g/h;
              f += e[j]*a[i][j];
            }
          hh=f/(h+h);
          for (j=0;j<=l;j++)
            { f=a[i][j];e[j]=g=e[j]-hh*f;
              for (k=0;k<=j;k++) a[j][k]-=(f*e[k]+g*a[i][k]);
            }
        }
      else e[i]=a[i][l];
      d[i]=h;
    }
                          // eigenvectors
#if (EVECTS == 1)
  d[0]=0.0;e[0]=0.0;
  for (i=0;i<n;i++)
    { l=i-1;
      if (d[i])
        { for (j=0;j<=l;j++)
            { g=0.0;
              for (k=0;k<=l;k++) g+=a[i][k]*a[k][j];
              for (k=0;k<=l;k++) a[k][j]-=g*a[k][i];
            }
        }
      d[i]=a[i][i];
      a[i][i]=1.0;
      for (j=0;j<=l;j++) a[j][i]=a[i][j]=0.0;
    }
#else 
  for (i=0;i<n;i++) d[i]=a[i][i];
#endif

  for (i=0;i<n-1;i++) e[i]=e[i+1];
  e[n-1]=0;
  // .......................................... complete transformation matrix
//  for (i=0;i<n;i++) for (j=i+1;j<n;j++) a[i][j]=a[j][i];
}

/*********************************************************************/
#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
int tqli(double *d, double *e, int n, double **z)
/***
 April 2005, Roger Melko, modified from Numerical Recipies in C v.2
 modified from www.df.unipi.it/~moruzzi/
 Diagonalizes a tridiagonal matrix: d[] is input as the diagonal elements,
 e[] as the off-diagonal.  If the eigenvalues of the tridiagonal matrix
 are wanted, input z as the identity matrix.  If the eigenvalues of the
 original matrix reduced by tred2 are desired, input z as the matrix
 output by tred2.  The kth column of z returns the normalized eigenvectors,
 corresponding to the eigenvalues output in d[k].
***/
{
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;

  for (l=0;l<n;l++) {
    iter=0;
    do { 
      for (m=l;m<n-1;m++) { 
	dd=fabs(d[m])+fabs(d[m+1]);
	if (fabs(e[m])+dd == dd) break;
      }
      if (m!=l) { 
	if (iter++ == 30) { 
	  cout <<"Too many iterations in tqli() \n";
	  return 0;
	}
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=sqrt((g*g)+1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) { 
	  f=s*e[i];
	  b=c*e[i];
	  if (fabs(f) >= fabs(g)) { 
	    c=g/f;r=sqrt((c*c)+1.0);
	    e[i+1]=f*r;
	    c *= (s=1.0/r);
	  }
	  else { 
	    s=f/g;r=sqrt((s*s)+1.0);
	    e[i+1]=g*r;
	    s *= (c=1.0/r);
	  }
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  p=s*r;
	  d[i+1]=g+p;
	  g=c*r-b;
#if (EVECTS == 1)
	  for (k=0;k<n;k++) { 
	    f=z[k][i+1];
	    z[k][i+1]=s*z[k][i]+c*f;
	    z[k][i]=c*z[k][i]-s*f;
	  }
#endif
	}
	d[l]=d[l]-p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m!=l);
  }
  return 1;
}
