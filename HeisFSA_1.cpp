// Elementary DMRG for Heisenberg S=1/2 chain, 
// infinite system algorithm to build chain
// ED with Lanczos
// finite system sweep
//
//Roger Melko May 9 2006

#include <blitz/array.h>
#include <fstream>

BZ_USING_NAMESPACE(blitz)

//template function prototypes
template<typename T> Array<T,2> reduceM2M2(const Array<T,4>&, const int);
void EigenValuesLAN(Array<double,4>&, Array<double,2>&, const int, double *);
void DMlargeEigen(Array<double,2>&, Array<double,2>&, const int, const int);

//block object
struct BLOCK {
  int size;    //# of sites
  Array<double,2> HAp;  //A' block hamiltonian
  Array<double,2> SzL;   //Sz(left) operator  
  Array<double,2> SmL;   //Sm(left) operator
  Array<double,2> SpL;   //Sp(left) operator
};

int main()
{
  int iter, NumI;
  int i1, i1p, i2, i2p, i3, i4; 
  int b1; 
  int m, st;     //# states
  int sites;     //# sites
  double Eval;   //Eigenvalue
  BLOCK blk;
  ofstream fout;
  char fname[7];
 
  fname[0] = '0'; fname[1] = '0'; fname[2] = '0'; 
  fname[3] = '.';
  fname[4] = 'E'; fname[5] = 'n'; fname[6] = 'v';

  cout<<"# states to keep: ";
  cin>>m;

  //Matrices
  Array<double,4> TSR(2,2,2,2);  //tensor product for Hab hamiltonian

  Array<double,4> Habcd(4,4,4,4); //superblock hamiltonian
  Array<double,2> HAB(4,4);        //new SYSTEM Hamiltonian
  Array<double,2> Psi(4,4); // ground state wavefunction
  Array<double,2> rhoTSR(4,4); // reduced density matrix
  Array<double,2> OO(m,4);   //the TRUNCATION matrix
  Array<double,2> OT(4,m);   // trasposed truncation matrix
  Array<double,2> Hl(4,m);   // the left half of new system H

  //create identity matrices
  Array<double,2> I2(2,2), I4(4,4), I2m(2*m,2*m);   
  I2=0.0; I4=0.0; I2m=0.0;
  for (b1=0; b1<2; b1++) I2(b1,b1)=1.0;
  for (b1=0; b1<4; b1++) I4(b1,b1)=1.0;
  for (b1=0; b1<(2*m); b1++) I2m(b1,b1)=1.0;
  Array<double,2> I2st(4,4);
  I2st = I4;

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
  Array<double,2> SzAB = reduceM2M2(TSR,2);

  TSR = Sm(i,k)*I2(j,l);
  Array<double,2> SmAB = reduceM2M2(TSR,2);

  TSR = Sp(i,k)*I2(j,l);
  Array<double,2> SpAB = reduceM2M2(TSR,2);

  HAB = H12;
  st = 2;     //start with a 2^2=4 state system
  sites = 2;

  /******infinite system algorithm loop first iteration**********/
  /******build system until number of desired states m **********/

  while(st <= m) {

   Habcd = HAB(i,k)*I2st(j,l) + I2st(i,k)*HAB(j,l) +
      SzAB(i,k)*SzAB(j,l)+ 0.5*SpAB(i,k)*SmAB(j,l) + 0.5*SmAB(i,k)*SpAB(j,l);
  
    EigenValuesLAN(Habcd,Psi,(4*st*st),&Eval);
     
    cout<<"sites: "<<2.0*sites;
    cout<<" Energy: "<<Eval/(2.0*sites)<<endl;
  
    rhoTSR = 0;
    for (i1=0; i1< 2*st; i1++)
      for (i1p=0; i1p< 2*st; i1p++)
        for (i2=0; i2< 2*st; i2++)
 	  rhoTSR(i1,i1p) += Psi(i1,i2)*Psi(i1p,i2); 

    OO.resize(2*st,2*st);
    OT.resize(2*st,2*st);
    DMlargeEigen(rhoTSR, OO, 2*st, 2*st);   
    for (i1=0; i1<2*st; i1++)
      for (i2=0; i2< 2*st; i2++)
        OT(i1,i2) = OO(i2,i1);  //Transpose

    blk.HAp.resize(2*st,2*st);
    //transform Operator Matrices to new basis
    Hl.resize(2*st,2*st);
    Hl = sum(HAB(i,k)*OT(k,j),k);   //Ha'
    blk.HAp = sum(OO(i,k)*Hl(k,j),k);   //(inner product)
    
    blk.SzL.resize(2*st,2*st);
    Hl = sum(SzAB(i,k)*OT(k,j),k);  
    blk.SzL = sum(OO(i,k)*Hl(k,j),k);
    
    blk.SpL.resize(2*st,2*st);
    Hl = sum(SpAB(i,k)*OT(k,j),k);  
    blk.SpL = sum(OO(i,k)*Hl(k,j),k);
    
    blk.SmL.resize(2*st,2*st);
    Hl = sum(SmAB(i,k)*OT(k,j),k);  
    blk.SmL = sum(OO(i,k)*Hl(k,j),k);

    // Add a single spin  (m*m*2*2 tensor)
    st *= 2;
    TSR.resize(st,2,st,2);
    
    TSR = blk.HAp(i,k)*I2(j,l) + blk.SzL(i,k)*Sz(j,l)+ 
      0.5*blk.SpL(i,k)*Sm(j,l) + 0.5*blk.SmL(i,k)*Sp(j,l) ;
    HAB.resize(2*st,2*st);            //Hamiltonian for next iteration
    HAB = reduceM2M2(TSR,st);
   //   cout<<HAB<<endl;

    SzAB.resize(2*st,2*st);  //Operators for next iteration
    TSR = I2st(i,k)*Sz(j,l);
    SzAB = reduceM2M2(TSR,st);

    SpAB.resize(2*st,2*st);
    TSR = I2st(i,k)*Sp(j,l);
    SpAB = reduceM2M2(TSR,st);

    SmAB.resize(2*st,2*st);
    TSR = I2st(i,k)*Sm(j,l);
    SmAB = reduceM2M2(TSR,st);  

    Habcd.resize(2*st,2*st,2*st,2*st);   //re-prepare superblock matrix
    Psi.resize(2*st,2*st);             //GS wavefunction
    rhoTSR.resize(2*st,2*st);
  
//    st *= 2;
    sites += 1;

    I2st.resize(4*st,4*st);    //redifine identity matrix
    I2st = 0.0;
    for (b1=0; b1<(4*st); b1++) I2st(b1,b1)=1.0;

  }//end while

  blk.size = sites-1;
  fname[2] = 48 + (blk.size)%10;          //some ASCII crap
  fname[1] = 48 + (blk.size-blk.size%10)/100;
  fname[0] = 48;
  fout.open(fname,ios::out);
  fout << blk.size <<endl;
  fout << blk.HAp ;
  fout << blk.SzL ;
  fout << blk.SpL ;
  fout << blk.SmL ;
  fout.close();
  
  OO.resize(m,(2*st));   
  OT.resize((2*st),m);  
  Hl.resize((2*st),m);   

  /*******************second+ iteration down here*************/

  for (iter = 1; iter<NumI; iter++){
    
    Habcd = HAB(i,k)*I2st(j,l) + I2st(i,k)*HAB(j,l) +
      SzAB(i,k)*SzAB(j,l)+ 0.5*SpAB(i,k)*SmAB(j,l) + 0.5*SmAB(i,k)*SpAB(j,l);

    EigenValuesLAN(Habcd,Psi,(4*st*st),&Eval);
    cout<<2*sites<<" "<<Eval/(2.0*sites)<<endl;
    
    rhoTSR = 0;
    for (i1=0; i1< 2*st; i1++)
      for (i1p=0; i1p< 2*st; i1p++)
        for (i2=0; i2< 2*st; i2++)
 	  rhoTSR(i1,i1p) += Psi(i1,i2)*Psi(i1p,i2); 

    DMlargeEigen(rhoTSR, OO, 2*st, m);   
    for (i1=0; i1<2*st; i1++)
      for (i2=0; i2< m; i2++)
        OT(i1,i2) = OO(i2,i1);  //Transpose

    st = m;

    if (iter == 1) blk.HAp.resize(m,m);
    //transform Operator Matrices to new basis
    Hl = sum(HAB(i,k)*OT(k,j),k);   //Ha'
    blk.HAp = sum(OO(i,k)*Hl(k,j),k);   //(inner product)
    
    if (iter == 1) blk.SzL.resize(m,m);
    Hl = sum(SzAB(i,k)*OT(k,j),k);  
    blk.SzL = sum(OO(i,k)*Hl(k,j),k);
    
    if (iter == 1) blk.SpL.resize(m,m);
    Hl = sum(SpAB(i,k)*OT(k,j),k);  
    blk.SpL = sum(OO(i,k)*Hl(k,j),k);
    
    if (iter == 1) blk.SmL.resize(m,m);
    Hl = sum(SmAB(i,k)*OT(k,j),k);  
    blk.SmL = sum(OO(i,k)*Hl(k,j),k);
    
   if (iter == 1) TSR.resize(m,2,m,2);
    // Add a single spin  (m*m*2*2 tensor)
    TSR = blk.HAp(i,k)*I2(j,l) + blk.SzL(i,k)*Sz(j,l)+ 
      0.5*blk.SpL(i,k)*Sm(j,l) + 0.5*blk.SmL(i,k)*Sp(j,l) ;
    if (iter == 1) HAB.resize(2*m,2*m);
    HAB = reduceM2M2(TSR,m);    

    if (iter == 1){
      SzAB.resize(2*m,2*m);  //Operators for next iteration
      TSR = I2m(i,k)*Sz(j,l);
      SzAB = reduceM2M2(TSR,m);
  
      SpAB.resize(2*m,2*m);
      TSR = I2m(i,k)*Sp(j,l);
      SpAB = reduceM2M2(TSR,m);
  
      SmAB.resize(2*m,2*m);
      TSR = I2m(i,k)*Sm(j,l);
      SmAB = reduceM2M2(TSR,m);
    }

   st = m;
   sites += 1;

   if (iter == 1){    //final resize
    I2st.resize(2*m,2*m);  
    I2st = I2m;

    Habcd.resize(2*m,2*m,2*m,2*m);   //re-prepare superblock matrix
    Psi.resize(2*m,2*m);             //GS wavefunction
    rhoTSR.resize(2*m,2*m);
    OO.resize(m,(2*m));
    OT.resize((2*m),m);
    Hl.resize((2*m),m);

   }

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

