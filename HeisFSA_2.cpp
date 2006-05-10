// Elementary DMRG for Heisenberg S=1/2 chain, 
// infinite system algorithm to build chain
// ED with Lanczos
// finite system sweep
//
//Roger Melko May 9 2006

#include <blitz/array.h>
#include <fstream>

BZ_USING_NAMESPACE(blitz)

//block object
struct BLOCK {
  int size;    //# of sites
  Array<double,2> HAp;  //A' block hamiltonian
  Array<double,2> SzL;   //Sz(left) operator  
  Array<double,2> SmL;   //Sm(left) operator
  Array<double,2> SpL;   //Sp(left) operator
};

//template function prototypes
template<typename T> Array<T,2> reduceM2M2(const Array<T,4>&, const int);
void EigenValuesLAN(Array<double,4>&, Array<double,2>&, const int, double *);
void DMlargeEigen(Array<double,2>&, Array<double,2>&, const int, const int);
void BlockWrite(BLOCK *, const int, char []);
void BlockRead(BLOCK *, const int, char []);

int main()
{
  int iter, NumI;
  int i1, i1p, i2, i2p, i3, i4; 
  int b1; 
  int m, st;     //# states
  int truncflag;
  int sites;     //# sites
  double Eval;   //Eigenvalue
  BLOCK blk;
  char fname[7];

  //initialize filename
  fname[0] = '.'; fname[1] = 48;  //ASCII for 0
  fname[4] = '.'; fname[5] = 'e'; 
  fname[6] = 'n'; //fname[7] = 'v';

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


  /******infinite system algorithm loop   ***********************/
  /******build system until number of desired states m **********/

  truncflag = 0;

  for (iter = 1; iter<NumI; iter++){

    Habcd = HAB(i,k)*I2st(j,l) + I2st(i,k)*HAB(j,l) +
      SzAB(i,k)*SzAB(j,l)+ 0.5*SpAB(i,k)*SmAB(j,l) + 0.5*SmAB(i,k)*SpAB(j,l);
  
    EigenValuesLAN(Habcd,Psi,(4*st*st),&Eval);
    
//     cout<<"sites: "<<2.0*sites;
//     if (truncflag == 0) cout<<" e ";
//     else cout<<" t ";
    cout<<2.0*sites<<" "<<1.0/(2.0*sites);
    cout<<" "<<Eval/(2.0*sites)<<endl;
    
    rhoTSR = 0;
    for (i1=0; i1< 2*st; i1++)
      for (i1p=0; i1p< 2*st; i1p++)
        for (i2=0; i2< 2*st; i2++)
 	  rhoTSR(i1,i1p) += Psi(i1,i2)*Psi(i1p,i2); 

    if (2*st <= m){     // NO TRUNCATION
      OO.resize(2*st,2*st);
      OT.resize(2*st,2*st);
      DMlargeEigen(rhoTSR, OO, 2*st, 2*st);   
      for (i1=0; i1<2*st; i1++)
	for (i2=0; i2< 2*st; i2++)
	  OT(i1,i2) = OO(i2,i1); 
      Hl.resize(2*st,2*st);
      blk.HAp.resize(2*st,2*st);
      blk.SzL.resize(2*st,2*st);
      blk.SpL.resize(2*st,2*st);
      blk.SmL.resize(2*st,2*st);
      st *= 2;
    }
    else {            // TRUNCATION
      if (truncflag == 0 || truncflag == 3){
	OO.resize(m,2*st);
	OT.resize(2*st,m);
	Hl.resize(2*st,m);
	truncflag ++; // 1 or 4
      }
      DMlargeEigen(rhoTSR, OO, 2*st, m);   
      for (i1=0; i1<2*st; i1++)
	for (i2=0; i2< m; i2++)
	  OT(i1,i2) = OO(i2,i1);
    }
    
    if (truncflag < 2){
      if (truncflag == 1) {
	truncflag = 2;	
	st = m;
	blk.HAp.resize(m,m);
	blk.SzL.resize(m,m);
	blk.SpL.resize(m,m);
	blk.SmL.resize(m,m);
      }
      
      TSR.resize(st,2,st,2);
    }

    //transform Operator Matrices to new basis
    Hl = sum(HAB(i,k)*OT(k,j),k);   //Ha'
    blk.HAp = sum(OO(i,k)*Hl(k,j),k);   //(inner product)
    
    Hl = sum(SzAB(i,k)*OT(k,j),k);  
    blk.SzL = sum(OO(i,k)*Hl(k,j),k);    
    
    Hl = sum(SpAB(i,k)*OT(k,j),k);  
    blk.SpL = sum(OO(i,k)*Hl(k,j),k);
    
    Hl = sum(SmAB(i,k)*OT(k,j),k);  
    blk.SmL = sum(OO(i,k)*Hl(k,j),k);
    
    TSR = blk.HAp(i,k)*I2(j,l) + blk.SzL(i,k)*Sz(j,l)+ 
      0.5*blk.SpL(i,k)*Sm(j,l) + 0.5*blk.SmL(i,k)*Sp(j,l) ;
    HAB.resize(2*st,2*st);            //Hamiltonian for next iteration
    HAB = reduceM2M2(TSR,st);
   //   cout<<HAB<<endl;
    
    if (truncflag < 3){
      if  (truncflag == 2) {
	truncflag = 3;
	I2st.resize(2*m,2*m);    //redifine identity matrix
	I2st = I2m;
      }
      else{
	I2st.resize(4*st,4*st);    //redifine I (still growing)
	I2st = 0.0;
	for (b1=0; b1<(4*st); b1++) I2st(b1,b1)=1.0;
      }
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
    }

    //fix this: inside BlockWrite function
    fname[3] = 48;// + (sites)%10;          //some ASCII crap
    fname[2] = 48;// + sites/10;
    BlockWrite(&blk,sites,fname);

    sites += 1;    

  }//end iteration


//   BlockRead(&blk,5,fname);

//   cout<<blk.SzL;

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

/*****************************************************************/
void BlockWrite(BLOCK *blk, const int sites, char fname[])
{
  ofstream fout;  

  blk->size = sites;

  fout.open(fname,ios::out);
  fout << blk->size <<endl;
  fout << blk->HAp ;
  fout << blk->SzL ;
  fout << blk->SpL ;
  fout << blk->SmL ;
  fout.close();

} //BlockWrite

/*****************************************************************/
void BlockRead(BLOCK *blk, const int sites, char fname[])
{
  ifstream fin;  
 
  blk->size = sites;
  fname[2] = 48 + (blk->size)%10;          //some ASCII crap
//  fname[1] = 48 + (blk->size-blk->size%10)/100;
//  fname[0] = 48;

  fin.open(fname,ios::in);
  fin >> blk->size; 
  fin >> blk->HAp ;
  fin >> blk->SzL ;
  fin >> blk->SpL ;
  fin >> blk->SmL ;
  fin.close();

//  cout<<blk->HAp;

}
