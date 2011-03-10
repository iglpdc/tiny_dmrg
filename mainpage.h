/**
 * @file mainpage.h
 *
 * @brief Just file with comments only to cheat Doxygen to make a nice
 * mainpage
 *
 * \mainpage
 *
 * @author Roger Melko 
 * @author Ivan Gonzalez
 * @date March 18th 2011
 *
 * @brief Elementry DMRG simulation for the Heisenberg chain; 
 *  \f$H= \sum_{i} (S^x_i S^x_{i+1} + S^y_i S^y_{i+1} +  S^z_i S^z_{i+1} ) \f$
 *
 * <ul> 
 *  <li> Begins with the "symmetric" infinite system algorithm to build up the chain
 *  <li> After this, a number of finite system algorithm sweeps are performed
 *  <li> The exact diagonalization performed with Lanczos
 *  <li> The output is the energy as a function of sweep
 *  <li> The code uses Blitz++ to handle tensors and matrices: see http://www.oonumerics.org/blitz/
 *  </ul>
 *
 * The outline of this document is the following:
 * <ul>
 * <li> \subpage start
 * <li> \subpage exercises   
 * <li> \subpage people   
 * </ul>
 *
 * \page start Getting Started
 *
 * \section reqs Requirements
 *
 * To be able to compile the code you just need any decent C++ compiler.
 * We have succesfully compiled the code with the following OS and
 * compilers:
 *
 * <ul> <li> Mac OS 10.4, 10.5, and 10.6 with the gcc compiler version >4.0 
 * (you just need to install the Developer Tools)
 * <li> Windows 7 with <a href="http://www.mingw.org/">MinGW</a> <li> Linux with with the gcc compiler
 * version >4.0 </ul>
 *
 * No libraries or other software are required, although having installed make and
 * some plotting software may be handy.
 *
 * \section download Download and Install
 *
 * You can get the source code in 
 * <a href="http://saskeram.cmt.uwaterloo.ca/dmrg_loo/dmrgloo.tar.gz">this link</a>. To uncompress the tarball do:
 *
 * \code $ tar -xzvf dmrgloo.tar.gz \endcode
 *
 * This will uncompress the source code in a directory called ./dmrgloo
 *
 * \section compile How to compile the code
 *
 * To compile the code just use in the directory that contains your source
 * files: 
 *
 * \code 
 * $ cd ./dmrgloo 
 * $ make 
 * \endcode
 *
 * If you do not have make installed, run this command instead:
 *
 * \code
 * $ g++ -c -O3 -I. tqli2.cpp tred3.cpp heisenberg.cpp densityMatrix.cpp lanczosDMRG.cpp
 * \endcode
 *
 * This will make a executable file called a.out that
 * implements the DMRG algorithm for the one-dimensional Heisenberg model.
 *
 * \section run Running the code
 *
 * To run the code do: 
 *
 * \code $ ./a.out \endcode
 *
 * \subsection params About the parameters
 *
 * This code is designed to be use as a didactical tool so it's far from
 * be optimized to do realistic DMRG calculations. While in a research
 * code one can easily keep hundreds of states in the truncated basis,
 * here we are restricted to keep a few tens, as the CPU resources used by
 * the code grow with the fourth power of the number of states kept.
 * Computational resources grow linearly with the other two parameters
 * (number of site of the chain and the number of sweeps). 
 *
 * In any decent laptop a choice of parameters could be:
 *
 * \code 
 * number of states to keep: 16 
 * number of sites in the chain: 16
 * number of FSA sweeps: 5 
 * \endcode
 *
 *
 * \page exercises Exercises
 * 
 * \section heisenberg Calculation of the energies of the Heisenberg model
 *
 * The first goal of the tutorial is to give you some background on the
 * implementation of the DMRG algorithm by looking at an actual code  
 * for the Heisenberg model.
 *
 * <ul> 
 * <li> Check out the documentation in this website and get familiar with the code.  
 * <li> Run the code and plot the energies as a function of the
 * site. <li> Compare the energy you obtain with the exact solution from
 * the Bethe Ansatz.  <li> Run the code with
 * different parameters and see the differences. </ul> 
 *
 * The exact result (Bethe Ansatz) for the infinite chain is \f$E_{0}=\frac{1}{4}-\ln 2\f$.
 *
 *
 *
 * \section entanglement Calculation of the entanglement entropy
 *
 * Implement a function to calculate the von Neumman entanglement
 * entropy at every DMRG step and print the result on screen as is done
 * for the energy. Remember that the von Neumann entanglement entropy is defined as:
 *
 * \f$S_{vN}(l)=-tr(\rho_{l}\ln\rho_{l})=-\sum_{i}\lambda_{i}\ln\lambda_{i}\f$,
 *
 * where \f$\lambda_{i}\f$ are the eigenvalues of the reduced density
 * matrix \f$\rho_{l}\f$.
 *
 * Plot the entanglement entropy and analyze your results. A detailed
 * explanation of the behaviour of the entanglement entropy is found in
 * <a href="http://prl.aps.org/abstract/PRL/v96/i10/e100603">Phys. Rev. Lett. 96,
 * 100603 (2006)</a>
 *
 *
 * \section ising Ising model in a tranverse field 
 *
 * Implement the DMRG algorithm for the one-dimensional S=1/2 Ising model
 * in a transverse magnetic field. The Hamiltonian of the Ising model in a
 * transverse field is:
 *
 * \f$H=J\sum_{i}S^{x}_{i}S^{x}_{i+1}+\Gamma\sum_{i}S^{z}_{i}\f$
 *
 * where \f$\vec{S}\f$ are the spin operators, \f$J<0\f$ is the
 * ferromagnetic coupling between neighboring spins, and \f$\Gamma\f$ is a
 * magnetic field. This model is important in condensed matter physics as
 * the simplest model having a quantum phase transition (QPT): at
 * \f$h_{c}=0.5\Gamma/J\f$ the ground state of the system changes from a
 * ferromagnetic to a paramagnetic phase.
 *
 * The main point of this exercise is to get used to the way the
 * hamiltonian is written within DMRG (splitting the whole chain in system and
 * environment). 
 *
 * <ul>
 * <li> Modify the heisenberg.cpp file, replacing the Hamiltonian by
 * the one of the Ising model.
 * <li> Think which operators you need to construct the Hamiltonian and
 * do the proper updates for them.
 * <li> Plot the energies in the center of the chain as a function to the magnetic field 
 * (you will need a new parameter corresponding to \f$h=\Gamma/J\f$).
 * <li> Plot the entanglement entropy as in the previous exercise. Do you
 * see something at the QPT?
 * </ul>
 *
 * \page people People
 *
 * Roger Melko, Ivan Gonzalez, Ann Kallin, and Kevin Resch
 *
 */
