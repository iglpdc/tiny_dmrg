/** 
 * @file heisenberg.cpp
 * @brief The main c++ file for the DMRG
 * 
 * @author Roger Melko 
 * @author Ivan Gonzalez
 * @date February March 7th, 2011
 *
 * @include mainpageForDoygen.cpp
 * 
 * @brief Elementry DMRG simulation for the Heisenberg chain; 
 *  \f$H= \sum_{ij} (S^x_i S^x_j + S^y_i S^y_j +  S^z_i S^z_j ) \f$
 *
 * <ul> 
 *  <li> Uses a "symmetric" infinite system algorithm to build chain 
 *  <li> Exact diagonalization performed with Lanczos
 *  <li> finite system sweep - symmetric build of L and R blocks
 *  </ul>
 */
#include "blitz/array.h"
#include "block.h"
#include "matrixManipulation.h"
#include "lanczosDMRG.h"
#include "densityMatrix.h"
#include "main_helpers.h"

int main()
{
    // Read some input 
    int numberOfHalfSweeps;
    int numberOfSites;    
    int m;
    std::cout<<"# states to keep: ";
    std::cin>>m;
    std::cout<<"System size : ";
    std::cin>>numberOfSites;
    std::cout<<"FSA sweeps : ";
    std::cin>>numberOfHalfSweeps;
    // done reading input

    /** \var BLOCK blkS \brief create the system block */
    BLOCK blkS;  //the system block
    ///create the environment block
    BLOCK blkE;  //the environment block
    //Matrices
    blitz::Array<double,4> TSR(2,2,2,2);  //tensor product for Hab hamiltonian

    blitz::Array<double,4> Habcd(4,4,4,4); //superblock hamiltonian
    blitz::Array<double,2> Psi(4,4);       // ground state wavefunction
    blitz::Array<double,2> rhoTSR(4,4);    // reduced density matrix
    blitz::Array<double,2> OO(m,4);        //the TRUNCATION matrix
    blitz::Array<double,2> OT(4,m);        // trasposed truncation matrix

    blitz::Array<double,2> HAp;  //A' block hamiltonian
    blitz::Array<double,2> SzB;   //Sz(left) operator  
    blitz::Array<double,2> SmB;   //Sm(left) operator
    blitz::Array<double,2> SpB;   //Sp(left) operator

    // create pauli matrices and identity matrix
    blitz::Array<double,2> Sz(2,2), Sp(2,2), Sm(2,2);
    Sz = 0.5, 0,
       0, -0.5;
    Sp = 0, 1.0,
       0, 0;
    Sm = 0, 0,
       1.0, 0;
    blitz::Array<double,2> I2=createIdentityMatrix(2);

    // tensor indices
    blitz::firstIndex i;    blitz::secondIndex j; 
    blitz::thirdIndex k;    blitz::fourthIndex l; 

    // build the Hamiltonian for two-sites only
    TSR = Sz(i,k)*Sz(j,l)+ 0.5*Sp(i,k)*Sm(j,l) + 0.5*Sm(i,k)*Sp(j,l);
    blkS.HAB.resize(4,4);
    blkS.HAB = reduceM2M2(TSR);

    TSR = Sz(i,k)*I2(j,l);
    blitz::Array<double,2> SzAB = reduceM2M2(TSR);

    TSR = Sm(i,k)*I2(j,l);
    blitz::Array<double,2> SmAB = reduceM2M2(TSR);

    TSR = Sp(i,k)*I2(j,l);
    blitz::Array<double,2> SpAB = reduceM2M2(TSR);
    // done building the Hamiltonian

    /**
     * Infinite system algorithm 
     */
    int st = 2;     //start with a 2^2=4 state system
    int sitesInSystem=2;     //# sites (SYSTEM)

    int truncflag = 0;
    int statesToKeepIFA=2; 
    blitz::Array<double,2> I2st=createIdentityMatrix(4);

    while (sitesInSystem <= (numberOfSites)/2 ) 
    {
	Habcd = blkS.HAB(i,k)*I2st(j,l)+ 
	    I2st(i,k)*blkS.HAB(j,l)+
	    SzAB(i,k)*SzAB(j,l)+0.5*SpAB(i,k)*SmAB(j,l)+0.5*SmAB(i,k)*SpAB(j,l);
      
	double groundStateEnergy=calculateGroundState(Habcd, Psi);

	printGroundStateEnergy(sitesInSystem, sitesInSystem, 
		groundStateEnergy);
	
	statesToKeepIFA= (2*statesToKeepIFA<=m)? 2*statesToKeepIFA : m;

        // decide whether you have to truncate or not 	
	if (2*st <= m){     // NO TRUNCATION

	    st *= 2;
	}
	else {            // TRUNCATION
	    if (truncflag == 0 || truncflag == 3){
		truncflag ++; // 1 or 4
	    }
	}
		
	// calculate the reduced density matrix and truncate 
	rhoTSR=calculateReducedDensityMatrix(Psi);

	OO.resize(statesToKeepIFA,rhoTSR.rows());
	OT.resize(rhoTSR.rows(),statesToKeepIFA);
	OO=truncateReducedDM(rhoTSR, statesToKeepIFA);   
	OT=OO.transpose(blitz::secondDim, blitz::firstDim);
	// end decision
	
	//transform the operators to new basis
	HAp.resize(statesToKeepIFA, statesToKeepIFA);
	SzB.resize(statesToKeepIFA, statesToKeepIFA);
	SpB.resize(statesToKeepIFA, statesToKeepIFA);
	SmB.resize(statesToKeepIFA, statesToKeepIFA);
	HAp=transformOperator(blkS.HAB, OT, OO);
	SzB=transformOperator(SzAB, OT, OO);
	SpB=transformOperator(SpAB, OT, OO);
	SmB=transformOperator(SmAB, OT, OO);
	
	// prepare Hamiltonian for next iter
	if (truncflag < 2){
	  if (truncflag == 1) {
	    truncflag = 2;	
	    st = m;
	  } // truncflag = 0
	}

	//Hamiltonian for next iteration
	TSR.resize(statesToKeepIFA,2,statesToKeepIFA,2);
	TSR = HAp(i,k)*I2(j,l) + SzB(i,k)*Sz(j,l)+ 
	  0.5*SpB(i,k)*Sm(j,l) + 0.5*SmB(i,k)*Sp(j,l) ;

	blkS.HAB.resize(2*st,2*st);            
	blkS.HAB = reduceM2M2(TSR);

	if (truncflag < 3){
	    int statesToKeep;

	    if  (truncflag == 2) {
		truncflag = 3;
		statesToKeep=2*m;
	    }
	    else{  // truncflag<2
		statesToKeep=4*st;
	    }

	    //redefine identity matrix
	    I2st.resize(statesToKeep, statesToKeep);    
	    I2st = createIdentityMatrix(statesToKeep);

	    //Operators for next iteration
	    SzAB.resize(2*st,2*st);  
	    TSR = I2st(i,k)*Sz(j,l);
	    SzAB = reduceM2M2(TSR);

	    SpAB.resize(2*st,2*st);
	    TSR = I2st(i,k)*Sp(j,l);
	    SpAB = reduceM2M2(TSR);

	    SmAB.resize(2*st,2*st);
	    TSR = I2st(i,k)*Sm(j,l);
	    SmAB = reduceM2M2(TSR);        

	    Habcd.resize(2*st,2*st,2*st,2*st);   //re-prepare superblock matrix
	    Psi.resize(2*st,2*st);             //GS wavefunction
	    rhoTSR.resize(2*st,2*st);
	} // if truncflag<3 

	sitesInSystem++;

	blkS.size = sitesInSystem;  //this is the size of the system block
	blkS.ISAwrite(sitesInSystem);
    
    }//end INFINITE SYSTEM ALGORITHM 
  
    std::cout<<"End of the infinite system algorithm\n";

    /**
     * Finite size algorithm 
     */
    {
	// find minimum size of the enviroment
	int minEnviromentSize=calculateMinEnviromentSize(m,numberOfSites);
     
	// start in the middle of the chain 
	int sitesInSystem = numberOfSites/2;
	blkS.FSAread(sitesInSystem,1);
      
	for (int halfSweep=0; halfSweep<numberOfHalfSweeps; halfSweep++)
	{
	    while (sitesInSystem <= numberOfSites-minEnviromentSize)
	    {
		int sitesInEnviroment = numberOfSites - sitesInSystem;

		// read the environment block from disk
		blkE.FSAread(sitesInEnviroment,halfSweep);

		// build the hamiltonian as a four-index tensor
		Habcd = blkE.HAB(i,k)*I2st(j,l)+I2st(i,k)*blkS.HAB(j,l)+
		SzAB(i,k)*SzAB(j,l)+
		0.5*SpAB(i,k)*SmAB(j,l)+0.5*SmAB(i,k)*SpAB(j,l);

		// calculate the ground state energy 
		double groundStateEnergy=calculateGroundState(Habcd, Psi);

		if (halfSweep%2 == 0) 
		    printGroundStateEnergy(sitesInSystem, sitesInEnviroment, 
			    groundStateEnergy);
		else 
		    printGroundStateEnergy(sitesInEnviroment, sitesInSystem, 
			    groundStateEnergy);
	    
		// calculate the reduced density matrix and truncate 
		rhoTSR=calculateReducedDensityMatrix(Psi);
	 
		blitz::Array<double,2> OO=truncateReducedDM(rhoTSR, m);   
		OT=OO.transpose(blitz::secondDim, blitz::firstDim);
	  
		// transform the operators to new basis
		HAp=transformOperator(blkS.HAB, OT, OO);
		SzB=transformOperator(SzAB, OT, OO);
		SpB=transformOperator(SpAB, OT, OO);
		SmB=transformOperator(SmAB, OT, OO);

		// add spin to the system block only
		TSR = HAp(i,k)*I2(j,l) + SzB(i,k)*Sz(j,l)+ 
		0.5*SpB(i,k)*Sm(j,l) + 0.5*SmB(i,k)*Sp(j,l);       
		blkS.HAB = reduceM2M2(TSR);
      
		sitesInSystem++;

		blkS.size = sitesInSystem;
		blkS.FSAwrite(sitesInSystem,halfSweep);
	    }// while

	    sitesInSystem = minEnviromentSize;
	    blkS.FSAread(sitesInSystem,halfSweep);

	}// for
    }  // end of the finite size algorithm
    return 0;
} // end main
