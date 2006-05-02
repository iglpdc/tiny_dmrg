// Elementary DMRG for Heisenberg S=1/2 chain, 
// infinite system algorithm
// ED with Lanczos
//
//Roger Melko Feb 26 2006

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

//template function prototypes
template<typename T> Array<T,2> reduceM2M2(const Array<T,4>&, const int);
void EigenValuesLAN(Array<double,4>&, Array<double,2>&, const int, double *);
void DMlargeEigen(Array<double,2>&, Array<double,2>&, const int, const int);

int main()
{
  int iter, NumI;
  int i1, i1p, i2, i2p, i3, i4; 
  int b1; 
  int m, st;      //# states
  double Eval;   //Eigenvalue

  cout<<"# states to keep: ";
  cin>>m;

  //Matrices
  Array<double,4> TSR(2,2,2,2);  //tensor product for Hab hamiltonian
  Array<double,4> TSR2m(2,m,2,m);

  Array<double,4> Habcd(4,4,4,4); //superblock hamiltonian
  Array<double,2> HAB(4,4);        //new SYSTEM Hamiltonian
  Array<double,2> Psi(4,4); // ground state wavefunction
  Array<double,2> rhoTSR(4,4); // reduced density matrix
  Array<double,2> OO(m,4);   //the TRUNCATION matrix
  Array<double,2> OT(4,m);   // trasposed truncation matrix
  Array<double,2> Hl(4,m);   // the left half of new system H

  Array<double,2> HAp(m,m);  //A' block hamiltonian
  Array<double,2> SzL(m,m);   //Sz(left) operator  
  Array<double,2> SmL(m,m);   //Sm(left) operator
  Array<double,2> SpL(m,m);   //Sp(left) operator

  //create identity matrices
  Array<double,2> I2(2,2), I4(4,4), Im(m,m), I2m(2*m,2*m);   
  I2=0.0; I4=0.0; I2m=0.0;
  for (b1=0; b1<2; b1++) I2(b1,b1)=1.0;
  for (b1=0; b1<4; b1++) I4(b1,b1)=1.0;
  for (b1=0; b1<m; b1++) Im(b1,b1)=1.0;
  for (b1=0; b1<(2*m); b1++) I2m(b1,b1)=1.0;

  // Create pauli matrices
  Array<double,2> Sz(2,2), Sp(2,2), Sm(2,2);
  Sz = 0.5, 0,
       0, -0.5;
  Sp = 0, 1.0,
       0, 0;
  Sm = 0, 0,
       1.0, 0;

  //tensor indices
  firstIndex i;    secondIndex j; 
  thirdIndex k;    fourthIndex l; 

  cout<<"Iterations : ";
  cin>>NumI;

  //create a tensor product: two-site Hamiltonian
  TSR = Sz(i,k)*Sz(j,l)+ 0.5*Sp(i,k)*Sm(j,l) + 0.5*Sm(i,k)*Sp(j,l) ;
  //write as 2D matrix in combined basis
  Array<double,2> H12 = reduceM2M2(TSR,2);
  //cout<<"H12 "<<H12<<endl;

  TSR = Sz(i,k)*I2(j,l);
  Array<double,2> SZab = reduceM2M2(TSR,2);
  
  TSR = Sm(i,k)*I2(j,l);
  Array<double,2> SMab = reduceM2M2(TSR,2);

  TSR = Sp(i,k)*I2(j,l);
  Array<double,2> SPab = reduceM2M2(TSR,2);

  HAB = H12;
  st = 4;     //start with a 2^2=4 state system

  /******infinite system algorithm loop first iteration**********/
  /******build system until number of desired states m **********/

  Array<double,2> SzAB(2*st,2*st); 
  Array<double,2> SpAB(2*st,2*st);
  Array<double,2> SmAB(2*st,2*st);

  while(st <= m) {

    Habcd = HAB(i,k)*I4(j,l) + I4(i,k)*HAB(j,l) + 
      SZab(i,k)*SZab(j,l)+ 0.5*SPab(i,k)*SMab(j,l) + 0.5*SMab(i,k)*SPab(j,l);
  
    EigenValuesLAN(Habcd,Psi,(st*st),&Eval);
     
    cout<<"Energy "<<Eval<<endl;
  
    HAp = HAB;   /*USE IF NO RG FOR 1st STEP*/
    SzL = SZab;
    SpL = SPab;
    SmL = SMab;

    // Add a single spin  (m*m*2*2 tensor)
    TSR.resize(st,2,st,2);
    TSR = HAp(i,k)*I2(j,l) + SzL(i,k)*Sz(j,l)+ 
      0.5*SpL(i,k)*Sm(j,l) + 0.5*SmL(i,k)*Sp(j,l) ;
    //cout<<TSR<<endl;
  

    HAB.resize(2*st,2*st);            //Hamiltonian for next iteration
    HAB = reduceM2M2(TSR,st);
   //   cout<<HAB<<endl;

    SzAB.resize(2*st,2*st);  //Operators for next iteration
    TSR = Im(i,k)*Sz(j,l);
    SzAB = reduceM2M2(TSR,st);

    SpAB.resize(2*st,2*st);
    TSR = Im(i,k)*Sp(j,l);
    SpAB = reduceM2M2(TSR,st);

    SmAB.resize(2*st,2*st);
    TSR = Im(i,k)*Sm(j,l);
    SmAB = reduceM2M2(TSR,st);  

    Habcd.resize(2*st,2*st,2*st,2*st);   //re-prepare superblock matrix
    Psi.resize(2*st,2*st);             //GS wavefunction
    rhoTSR.resize(2*st,2*st);
    OO.resize(st,(2*st));   
    OT.resize((2*st),st);  
    Hl.resize((2*st),st);   
  
    st *= 2;

  }//end while

  /*******************second+ iteration down here*************/

  for (iter = 1; iter<NumI; iter++){
    
    Habcd = HAB(i,k)*I2m(j,l) + I2m(i,k)*HAB(j,l) +
      SzAB(i,k)*SzAB(j,l)+ 0.5*SpAB(i,k)*SmAB(j,l) + 0.5*SmAB(i,k)*SpAB(j,l);

    EigenValuesLAN(Habcd,Psi,(2*m)*(2*m),&Eval);
//    cout<<iter<<" Energy "<<Eval<<" "<<Eval/(2*iter+4)<<endl;
    cout<<1.0/(2.0*iter+4.0)<<" "<<Eval/(2.0*iter+4.0)<<endl;
    
    rhoTSR = 0;
    for (i1=0; i1<2*m; i1++)
      for (i1p=0; i1p<2*m; i1p++)
        for (i2=0; i2<2*m; i2++)
 	  rhoTSR(i1,i1p) += Psi(i1,i2)*Psi(i1p,i2); 


    DMlargeEigen(rhoTSR, OO, (2*m), m);   
    for (i1=0; i1<2*m; i1++)
      for (i2=0; i2<m; i2++)
        OT(i1,i2) = OO(i2,i1);  //Transpose

    //transform Operator Matrices to new basis
    Hl = sum(HAB(i,k)*OT(k,j),k);   //Ha'
    HAp = sum(OO(i,k)*Hl(k,j),k);   //(inner product)
    
    Hl = sum(SzAB(i,k)*OT(k,j),k);  
    SzL = sum(OO(i,k)*Hl(k,j),k);
    
    Hl = sum(SpAB(i,k)*OT(k,j),k);  
    SpL = sum(OO(i,k)*Hl(k,j),k);
    
    Hl = sum(SmAB(i,k)*OT(k,j),k);  
    SmL = sum(OO(i,k)*Hl(k,j),k);
    
    // Add a single spin  (m*m*2*2 tensor)
    TSR = HAp(i,k)*I2(j,l) + SzL(i,k)*Sz(j,l)+ 
      0.5*SpL(i,k)*Sm(j,l) + 0.5*SmL(i,k)*Sp(j,l) ;
    HAB = reduceM2M2(TSR,m);    

  }//iter

  return 0;
}

/*****************************************************************/
template<typename T>
Array<T,2> reduceM2M2(const Array<T,4>& T2, const int m)
/****reduces a m,2,m,2 matrix to a 2m*2m matrix
     over the direct product basis of 2 spins*************/
{
  int r,c;
  Array<double,2> M4(2*m,2*m);

  r=0;
  for (int a1=0; a1<m; a1++)
    for (int a2=0; a2<2; a2++){
      c=0;
      for (int a3=0; a3<m; a3++)
	for (int a4=0; a4<2; a4++){
	  M4(r,c) = T2(a1,a2,a3,a4);
	  c++;
	}
      r++;    
    }
  
  return M4;

}

