/** 
 * @file heisenberg.cpp
 * @brief The main c++ file for the DMRG
 *
 * 
 * @author Roger Melko 
 * @author Ivan Gonzalez
 * @date $Date$
 *
 * $Revision$ 
 *
 * Elementary DMRG simulation for the Heisenberg chain; 
 *  \f$H= \sum_{ij} (S^x_i S^x_j + S^y_i S^y_j +  S^z_i S^z_j ) \f$
 *
 * <ul> 
 *  <li> Begins with the "symmetric" infinite system algorithm to build 
 *  up the chain
 *  <li> After this, a number of finite system algorithm sweeps are performed
 *  <li> The exact diagonalization performed with Lanczos
 *  <li> The output is the energy as a function of sweep
 *  <li> The code uses Blitz++ to handle tensors and matrices: see http://www.oonumerics.org/blitz/
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
    // Read some input from user
    int numberOfHalfSweeps;
    int numberOfSites;    
    int m;
    std::cout<<"Enter the number states to keep: ";
    std::cin>>m;
    std::cout<<"Enter the number of sites in the chain: ";
    std::cin>>numberOfSites;
    std::cout<<"Enter the number of FSA sweeps : ";
    std::cin>>numberOfHalfSweeps;

    Block system;   //create the system block
    Block environ;  //create the environment block

    //Below we declare the Blitz++ matrices used by the program
    blitz::Array<double,4> TSR(2,2,2,2);  //tensor product for Hab hamiltonian

    blitz::Array<double,4> Habcd(4,4,4,4); // superblock hamiltonian
    blitz::Array<double,2> Psi(4,4);       // ground state wavefunction
    blitz::Array<double,2> reducedDM(4,4); // reduced density matrix
    blitz::Array<double,2> OO(m,4);        // the truncation matrix
    blitz::Array<double,2> OT(4,m);        // transposed truncation matrix

    blitz::Array<double,2> HAp;   //A' block hamiltonian
    blitz::Array<double,2> SzB;   //Sz(left) operator  
    blitz::Array<double,2> SmB;   //Sm(left) operator
    blitz::Array<double,2> SpB;   //Sp(left) operator

    // create the pauli matrices and the 2x2 identity matrix
    blitz::Array<double,2> sigma_z(2,2), sigma_p(2,2), sigma_m(2,2);
    sigma_z = 0.5, 0,
         0, -0.5;
    sigma_p = 0, 1.0,
         0, 0;
    sigma_m = 0, 0,
         1.0, 0;
    blitz::Array<double,2> I2=createIdentityMatrix(2);

    // declare tensor indices according to Blitz++ convention
    blitz::firstIndex i;    blitz::secondIndex j; 
    blitz::thirdIndex k;    blitz::fourthIndex l; 

    // build the Hamiltonian for two-sites only
    TSR = sigma_z(i,k)*sigma_z(j,l)+ 0.5*sigma_p(i,k)*sigma_m(j,l) + 
	0.5*sigma_m(i,k)*sigma_p(j,l);
    system.blockH.resize(4,4);
    system.blockH = reduceM2M2(TSR);

    TSR = sigma_z(i,k)*I2(j,l);
    blitz::Array<double,2> SzAB = reduceM2M2(TSR);

    TSR = sigma_m(i,k)*I2(j,l);
    blitz::Array<double,2> SmAB = reduceM2M2(TSR);

    TSR = sigma_p(i,k)*I2(j,l);
    blitz::Array<double,2> SpAB = reduceM2M2(TSR);
    // done building the Hamiltonian

    /**
     * Infinite system algorithm: build the Hamiltonian from 2 to N-sites
     */
    int statesToKeepIFA=2;   //start with a 2^2=4 state system
    int sitesInSystem=2;     //# sites (SYSTEM)

    blitz::Array<double,2> I2st=createIdentityMatrix(4);

    while (sitesInSystem <= (numberOfSites)/2 ) 
    {
        Habcd = system.blockH(i,k)*I2st(j,l)+ 
            I2st(i,k)*system.blockH(j,l)+
            SzAB(i,k)*SzAB(j,l)+0.5*SpAB(i,k)*SmAB(j,l)+0.5*SmAB(i,k)*SpAB(j,l);

        double groundStateEnergy=calculateGroundState(Habcd, Psi);

        printGroundStateEnergy(sitesInSystem, sitesInSystem, groundStateEnergy);

        statesToKeepIFA= (2*statesToKeepIFA<=m)? 2*statesToKeepIFA : m;

        // calculate the reduced density matrix and truncate 
        reducedDM=calculateReducedDensityMatrix(Psi);

        OO.resize(statesToKeepIFA,reducedDM.rows()); //resize transf. matrix
        OT.resize(reducedDM.rows(),statesToKeepIFA); // and its inverse
        OO=truncateReducedDM(reducedDM, statesToKeepIFA); //get transf. matrix 
        OT=OO.transpose(blitz::secondDim, blitz::firstDim); //and its inverse

        //transform the operators to new basis
        HAp.resize(statesToKeepIFA, statesToKeepIFA);
        SzB.resize(statesToKeepIFA, statesToKeepIFA);
        SpB.resize(statesToKeepIFA, statesToKeepIFA);
        SmB.resize(statesToKeepIFA, statesToKeepIFA);
        HAp=transformOperator(system.blockH, OT, OO);
        SzB=transformOperator(SzAB, OT, OO);
        SpB=transformOperator(SpAB, OT, OO);
        SmB=transformOperator(SmAB, OT, OO);

        //Hamiltonian for next iteration
        TSR.resize(statesToKeepIFA,2,statesToKeepIFA,2);
        TSR = HAp(i,k)*I2(j,l) + SzB(i,k)*sigma_z(j,l)+ 
            0.5*SpB(i,k)*sigma_m(j,l) + 0.5*SmB(i,k)*sigma_p(j,l) ;

        system.blockH.resize(2*statesToKeepIFA,2*statesToKeepIFA);            
        system.blockH = reduceM2M2(TSR);

	int statesToKeep= (2*statesToKeepIFA<=m)? 4*statesToKeepIFA : 2*m;

	//redefine identity matrix
	I2st.resize(statesToKeep, statesToKeep);    
	I2st = createIdentityMatrix(statesToKeep);

	//Operators for next iteration
	SzAB.resize(2*statesToKeepIFA,2*statesToKeepIFA);  
	TSR = I2st(i,k)*sigma_z(j,l);
	SzAB = reduceM2M2(TSR);

	SpAB.resize(2*statesToKeepIFA,2*statesToKeepIFA);
	TSR = I2st(i,k)*sigma_p(j,l);
	SpAB = reduceM2M2(TSR);

	SmAB.resize(2*statesToKeepIFA,2*statesToKeepIFA);
	TSR = I2st(i,k)*sigma_m(j,l);
	SmAB = reduceM2M2(TSR);        

	Habcd.resize(2*statesToKeepIFA,2*statesToKeepIFA,2*statesToKeepIFA,2*statesToKeepIFA);   //re-prepare superblock matrix
	Psi.resize(2*statesToKeepIFA,2*statesToKeepIFA);             //GS wavefunction
	reducedDM.resize(2*statesToKeepIFA,2*statesToKeepIFA);
	
        sitesInSystem++;

        system.size = sitesInSystem;  //this is the size of the system block
        system.ISAwrite(sitesInSystem);

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
        system.FSAread(sitesInSystem,1);

        for (int halfSweep=0; halfSweep<numberOfHalfSweeps; halfSweep++)
        {
            while (sitesInSystem <= numberOfSites-minEnviromentSize)
            {
                int sitesInEnviroment = numberOfSites - sitesInSystem;

                // read the environment block from disk
                environ.FSAread(sitesInEnviroment,halfSweep);

                // build the hamiltonian as a four-index tensor
                Habcd = environ.blockH(i,k)*I2st(j,l)+
		    I2st(i,k)*system.blockH(j,l)+
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
                reducedDM=calculateReducedDensityMatrix(Psi);

                blitz::Array<double,2> OO=truncateReducedDM(reducedDM, m);   
                OT=OO.transpose(blitz::secondDim, blitz::firstDim);

                // transform the operators to new basis
                HAp=transformOperator(system.blockH, OT, OO);
                SzB=transformOperator(SzAB, OT, OO);
                SpB=transformOperator(SpAB, OT, OO);
                SmB=transformOperator(SmAB, OT, OO);

                // add spin to the system block only
                TSR = HAp(i,k)*I2(j,l) + SzB(i,k)*sigma_z(j,l)+ 
                    0.5*SpB(i,k)*sigma_m(j,l) + 0.5*SmB(i,k)*sigma_p(j,l);       
                system.blockH = reduceM2M2(TSR);

                sitesInSystem++;

                system.size = sitesInSystem;
                system.FSAwrite(sitesInSystem,halfSweep);
            }// while

            sitesInSystem = minEnviromentSize;
            system.FSAread(sitesInSystem,halfSweep);

        }// for
    }  // end of the finite size algorithm
    return 0;
} // end main
