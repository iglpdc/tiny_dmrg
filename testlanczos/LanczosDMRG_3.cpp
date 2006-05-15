// Lanczos routine for performing ED in 
// DMRG
//
//Roger Melko Feb. 23 2006

#define STARTIT 4  //iteration which diagonz. begins

#include <cmath>
#include <iomanip>
#include <blitz/array.h>
#include <fstream>

BZ_USING_NAMESPACE(blitz)

//template function prototypes
int LanczosED(Array<double,2>&, const int);
void Normalize(Array<double,1>&, const int);
int tqli2(Array<double,1>&, Array<double,1>&, int, Array<double,2>&, const int);
void EigenValues(double **, double *, const int );


int main()
{
  int i,j, dim; //linear dimension of Ham
  int min;
  int lrt;
  double Eex;
 
  //Matrices
  Array<double,2> Ham2(1,1);
  
  cout<<"H lin dim: ";
  cin>>dim;

  Ham2.resize(dim,dim);

  //debug read in matrix
  ifstream fin;
  ofstream fout;

  fin.open("Ham.dat",ios::in);
  fin >> Ham2 ;
  fin.close();


//  for (i=0; i<dim; i++)
//    for (j=0; j<dim; j++){
//      Ham2(i,j) = 0.1*(rand()%100);
//      Eex = rand()%2;
//      if (Eex == 0) Ham2(i,j) *= -1.0000001;
//    }
//  for (i=0; i<dim; i++)
//   for (j=i+1; j<dim; j++)
//     Ham2(i,j) = Ham2(j,i); 
//   cout<<Ham2<<endl;

  fout.open("Ham.dat",ios::out);
  fout << Ham2;
  fout.close();
  

  //householderdiagonalize tri-di matrix
  double *dd = new double[dim];
  double **Hmatr = new double*[dim];  //dynamic 2D array
  for (i=0; i < dim; i++)           //for diagonalization
    Hmatr[i] = new double[dim];

  for (i=0; i<dim; i++)
    for (j=0; j<dim; j++)  
      Hmatr[i][j] = Ham2(i,j);
  
  EigenValues(Hmatr,dd,dim);
  
  min = 0;
  for (i=1;i< dim ;i++)
    if (dd[i] < dd[min])  min = i;
  Eex = dd[min];

 cout<<"EigenVec : ";
 for (i=0;i< dim; i++)
   cout<<setprecision(12)<<Hmatr[i][min]<<", ";
 cout<<endl;  

//   cout<<Ham2<<endl;

  for (int i=0; i < dim; i++)  //free up your memory
    delete[] Hmatr[i];
  delete [] Hmatr;
  delete [] dd;

  cout<<"AA \n";
  lrt = LanczosED(Ham2, dim);
  cout<<"BB \n";


  cout<<setprecision(12)<<"Exact e0 "<<Eex<<"\n";

  return 0;
	
}

/******************************************************************/
int LanczosED(Array<double,2>& Ham, const int N)
/*******************************************************************
 ** Reduces the Hamiltonian Matrix to a tri-diagonal form  ********/
{
  int ii, jj;
  int iter, MAXiter, EViter;
  int min;
  int Lexit;
  double Norm;
  double E0;

  int LIT;   //max number of Lanczos iterations
  LIT = 100;
  
  //Matrices
  Array<double,1> V0(N);  
  Array<double,1> Vorig(N);
  Array<double,1> V1(N);  //Ground state vector
  Array<double,1> V2(N);
  Array<double,1> Psi(N);  //return eigenvector
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
  
  iter = 0;
  //initialize 
  for (ii=0; ii<N; ii++){
    Vorig(ii) = rand()%10*0.1;  //random starting vec
    Norm = rand()%2;
    if (Norm == 0) Vorig(ii) *= -1.0000001;
  }
  Normalize(Vorig,N);  

  for (EViter = 0; EViter < 2; EViter++) {//0=get E0 converge, 1=get eigenvec

    iter = 0;
    V0 = Vorig;
    
    if (EViter == 1) Psi = V0*(Hmatrix(0,min));
    
    V1 = 0;
    beta(0)=0;  //beta_0 not defined
    
    V1 = sum(Ham(i,j)*V0(j),j); // V1 = H |V0> 
    
    alpha(0) = 0;
    for (ii=0; ii<N; ii++) alpha(0) += V0(ii)*V1(ii);
    
    V1 = V1- alpha(0)*V0;
    Norm = 0;
    for (ii=0; ii<N; ii++) Norm += V1(ii)*V1(ii);
    beta(1) = sqrt(Norm);

    V1 = V1/beta(1);

    if (EViter == 1) Psi += V1*(Hmatrix(1,min));
    
    // done 0th iteration
    
    Lexit = 0;   //exit flag
    E0 = 1.0;    //previous iteration GS eigenvalue
    while(Lexit != 1){
      
      iter++;
      
      V2 = sum(Ham(i,j)*V1(j),j); // V2 = H |V1>
      //V2 -= beta(iter)*V0;
      
      alpha(iter) = 0;
      for (ii=0; ii<N; ii++) alpha(iter) += V1(ii)*V2(ii);
      
      V2 = V2- alpha(iter)*V1 -  beta(iter)*V0;
      Norm = 0;
      for (ii=0; ii<N; ii++) Norm += V2(ii)*V2(ii);
      beta(iter+1) = sqrt(Norm);
      
      V2 = V2/beta(iter+1);

      if (EViter == 1) {Psi += V2*(Hmatrix(iter+1,min));
	//cout<<Psi<<" S \n";
      }
      
      V0 = V1;
      V1 = V2;
      
      if (iter > STARTIT && EViter == 0){
	
	//diagonalize tri-di matrix
	d(0) = alpha(0);
	for (ii=1;ii<=iter;ii++){
	  d(ii) = alpha(ii);
	  e(ii-1) = beta(ii);
	}
	e(iter) = 0;
	
	nn = iter+1;
	rtn = tqli2(d,e,nn,Hmatrix,0);
	
	min = 0;
	for (ii=1;ii<=iter;ii++)
	  if (d(ii) < d(min))  min = ii;
	
	if ( (E0 - d(min)) < 1E-8) {
	  Lexit = 1;
	  cout<<"Lanc :"<<iter<<" ";
	  cout<<setprecision(12)<<d(min)<<"\n";     
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
      for (ii=1;ii<=iter;ii++){
	d(ii) = alpha(ii);
	e(ii-1) = beta(ii);
      }
      e(iter) = 0;
      //calculate eigenvector
      Hmatrix = 0;
      for (ii=0;ii<=iter;ii++)
	Hmatrix(ii,ii) = 1.0; //identity matrix
      nn = iter+1;
      rtn = tqli2(d,e,nn,Hmatrix,1);
      min = 0;
      for (ii=1;ii<=iter;ii++)
	if (d(ii) < d(min))  min = ii;
    }
    
  }//repeat (EViter) to transfrom eigenvalues H basis
  
  Normalize(Psi,N);

  cout<<Psi<<" Psi \n";

  double Etemp;
  Etemp = 0;
  V2 = sum(Ham(i,j)*Psi(j),j);
  for (ii=0;ii<N;ii++)
    Etemp += Psi(ii)*V2(ii);
//    cout<<ii<<" "<<V2(ii)/Psi(ii)<<" EVdiv \n";
  cout<<Etemp<<" EVdiv \n";

  return 0;
  
} //LanczosED

/*********************************************************************/
void Normalize(Array<double,1>& V, const int N) 
/*********************************************************************
    Normalize the input vector (length N)
***/
{
  double norm;
  int i;

  norm = 0.0;             
  for (i=0; i<N; i++)
    norm += V(i)*V(i); //  <V|V>
  norm = sqrt(norm);

  for (i=0; i<N; i++)
    V(i) /= norm;
}


/*********************************************************************/
#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
int tqli2(Array<double,1>& d, Array<double,1>& e, int n, Array<double,2>& z, const int Evects)
/***
 April 2005, Roger Melko, modified from Numerical Recipies in C v.2
 modified from www.df.unipi.it/~moruzzi/
 Diagonalizes a tridiagonal matrix: d[] is input as the diagonal elements,
 e[] as the off-diagonal.  If the eigenvalues of the tridiagonal matrix
 are wanted, input z as the identity matrix.  If the eigenvalues of the
 original matrix reduced by tred2 are desired, input z as the matrix
 output by tred2.  The kth column of z returns the normalized eigenvectors,
 corresponding to the eigenvalues output in d[k].
 Feb 23 2005: modified to use Blitz++ arrays
***/
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
