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
    Block env;  //create the environment block

    //Below we declare the Blitz++ matrices used by the program
    blitz::Array<double,4> TSR(2,2,2,2);   //tensor product for Hab hamiltonian

    blitz::Array<double,4> Habcd(4,4,4,4); // superblock hamiltonian
    blitz::Array<double,2> Psi(4,4);       // ground state wavefunction
    blitz::Array<double,2> reducedDM(4,4); // reduced density matrix
    blitz::Array<double,2> OO(m,4);        // the truncation matrix
    blitz::Array<double,2> OT(4,m);        // transposed truncation matrix

    blitz::Array<double,2> blockH_p;       //block hamiltonian after transform.
    blitz::Array<double,2> S_z_p;          //S_z operator after transformation  
    blitz::Array<double,2> S_m_p;          //S_m operator after transformation
    blitz::Array<double,2> S_p_p;          //S_p operator after transformation

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
    blitz::Array<double,2> S_z = reduceM2M2(TSR);

    TSR = sigma_m(i,k)*I2(j,l);
    blitz::Array<double,2> S_m = reduceM2M2(TSR);

    TSR = sigma_p(i,k)*I2(j,l);
    blitz::Array<double,2> S_p = reduceM2M2(TSR);
    // done building the Hamiltonian

    /**
     * Infinite system algorithm: build the Hamiltonian from 2 to N-sites
     */
    int statesToKeep=2;      //start with a 2^2=4 state system
    int sitesInSystem=2;     //# sites in the system block

    blitz::Array<double,2> I2st=createIdentityMatrix(4);

    while (sitesInSystem <= (numberOfSites)/2 ) 
    {
	// build the hamiltonian as a four-index tensor
        Habcd = system.blockH(i,k)*I2st(j,l)+ 
            I2st(i,k)*system.blockH(j,l)+
            S_z(i,k)*S_z(j,l)+0.5*S_p(i,k)*S_m(j,l)+0.5*S_m(i,k)*S_p(j,l);

	// calculate the ground state energy 
        double groundStateEnergy=calculateGroundState(Habcd, Psi);

        printGroundStateEnergy(sitesInSystem, sitesInSystem, groundStateEnergy);

	// increase the number of states if you are not at m yet
        statesToKeep= (2*statesToKeep<=m)? 2*statesToKeep : m;

        // calculate the reduced density matrix and truncate 
        reducedDM=calculateReducedDensityMatrix(Psi);

        OO.resize(statesToKeep,reducedDM.rows()); //resize transf. matrix
        OT.resize(reducedDM.rows(),statesToKeep); // and its inverse
        OO=truncateReducedDM(reducedDM, statesToKeep); //get transf. matrix 
        OT=OO.transpose(blitz::secondDim, blitz::firstDim); //and its inverse

        //transform the operators to new basis
        blockH_p.resize(statesToKeep, statesToKeep);
        S_z_p.resize(statesToKeep, statesToKeep);
        S_p_p.resize(statesToKeep, statesToKeep);
        S_m_p.resize(statesToKeep, statesToKeep);
        blockH_p=transformOperator(system.blockH, OT, OO);
        S_z_p=transformOperator(S_z, OT, OO);
        S_p_p=transformOperator(S_p, OT, OO);
        S_m_p=transformOperator(S_m, OT, OO);

        //Hamiltonian for next iteration
        TSR.resize(statesToKeep,2,statesToKeep,2);
        TSR = blockH_p(i,k)*I2(j,l) + S_z_p(i,k)*sigma_z(j,l)+ 
            0.5*S_p_p(i,k)*sigma_m(j,l) + 0.5*S_m_p(i,k)*sigma_p(j,l) ;

        system.blockH.resize(2*statesToKeep,2*statesToKeep);            
        system.blockH = reduceM2M2(TSR);

	//redefine identity matrix
	int statesToKeepNext= (2*statesToKeep<=m)? 4*statesToKeep : 2*m;
	I2st.resize(statesToKeepNext, statesToKeepNext);    
	I2st = createIdentityMatrix(statesToKeepNext);

	//redefine the operators for next iteration
	S_z.resize(2*statesToKeep,2*statesToKeep);  
	TSR = I2st(i,k)*sigma_z(j,l);
	S_z = reduceM2M2(TSR);

	S_p.resize(2*statesToKeep,2*statesToKeep);
	TSR = I2st(i,k)*sigma_p(j,l);
	S_p = reduceM2M2(TSR);

	S_m.resize(2*statesToKeep,2*statesToKeep);
	TSR = I2st(i,k)*sigma_m(j,l);
	S_m = reduceM2M2(TSR);        

	// re-prepare superblock matrix, wavefunction and reduced DM
	Habcd.resize(2*statesToKeep,2*statesToKeep,2*statesToKeep,2*statesToKeep);   
	Psi.resize(2*statesToKeep,2*statesToKeep);             
	reducedDM.resize(2*statesToKeep,2*statesToKeep);

	// make the system one site larger and save it
        system.size = ++sitesInSystem;  
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
                env.FSAread(sitesInEnviroment,halfSweep);

                // build the hamiltonian as a four-index tensor
                Habcd = env.blockH(i,k)*I2st(j,l)+
		    I2st(i,k)*system.blockH(j,l)+
                    S_z(i,k)*S_z(j,l)+
                    0.5*S_p(i,k)*S_m(j,l)+0.5*S_m(i,k)*S_p(j,l);

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
                blockH_p=transformOperator(system.blockH, OT, OO);
                S_z_p=transformOperator(S_z, OT, OO);
                S_p_p=transformOperator(S_p, OT, OO);
                S_m_p=transformOperator(S_m, OT, OO);

                // add spin to the system block only
                TSR = blockH_p(i,k)*I2(j,l) + S_z_p(i,k)*sigma_z(j,l)+ 
                    0.5*S_p_p(i,k)*sigma_m(j,l) + 0.5*S_m_p(i,k)*sigma_p(j,l);       
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
