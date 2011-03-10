/**
 * @file mainpage.h
 *
 * @brief Just file with only comments to cheat Doxygen to make a nice
 * mainpage
 *
 * \mainpage
 *
 * @author Roger Melko 
 * @author Ivan Gonzalez
 * @date March 18th 2011
 *
 * @brief Elementry DMRG simulation for the Heisenberg chain; 
 *  \f$H= \sum_{ij} (S^x_i S^x_j + S^y_i S^y_j +  S^z_i S^z_j ) \f$
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
 * <ul> <li> Mac OS 10.4, 10.5, and 10.6 with the gcc compiler version
 * >4.0 <li> Windows 7 with MinGW <li> Linux with with the gcc compiler
 * version >4.0 </ul>
 *
 * No libraries or other software are required. Although having make and
 * some plotting software may be handy.
 *
 * \section download Download and install
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
 * \code $ cd ./dmrgloo $ make \endcode
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
 * The first goal of the tutorial is to give you some background on the
 * implementation of the DMRG algorithm. 
 *
 * <ul> <li> Check out the different pages in this website and get
 * familiar with the code.  <li> Plot the energies as a function of the
 * site <li> Compare the energy you obtained with the exact solution from
 * the Bethe Ansatz (put the number here).  <li> Run the code with
 * different parameters and see the differences </ul> 
 *
 * After you have played with the code a bit, we propose to implement to
 * new features on it. 
 * 
 * \section entanglement Calculation of the entanglement entropy
 *
 * Implement a function to calculate the (von Neumman) entanglement
 * entropy at every DMRG step and print the result on screen as is done
 * with the energy.
 *
 * Remember that the (von Neumann) entanglement entropy is defined as \f$
 * S_{vN}(l)=-tr(\rho_{l}\ln\rho_{l})=-\sum_{i}\lambda_{i}\ln\lambda_{i}\f$
 * where \f$\lambda_{i}\f$ are the eigenvalues of the reduced density
 * matrix \f$\rho_{l}\f$.
 *
 * Plot the entanglement entropy and analyze your results. A detailed
 * explanation of the behaviour of the entanglement entropy in this paper
 * (cite Affleck's paper)
 *
 * \section ising Ising model in a tranverse field 
 *
 * Implement the DMRG algorithm for the one-dimensional S=1/2 Ising model
 * in a tranverse field.  The Hamiltonian of the Ising model in a
 * transverse field is:
 *
 * \f$H=J\sum_{i}S^{x}_{i}S^{x}_{i+1}+\Gamma\sum_{i}S^{z}_{i}\f$
 *
 * where \f$\vec{S}\f$ are the spin operators, \f$J<0\f$ is the
 * ferromagnetic coupling between neighboring spins, and \f$\Gamma\f$ is a
 * magnetic field. This model is important in condensed matter physics as
 * it's one the simplest models having a quantum phase transition (QPT)): at
 * \f$h_{c}=0.5\Gamma/J\f$ the ground state of the system changes from a
 * ferromagnetic to a paramagnetic phase.
 *
 * The main point of this exercise is to get used with the way the
 * hamiltonian is written in DMRG (splitting the whole chain in system and
 * environment). 
 *
 * <ul>
 * <li> Modify the heisenberg.cpp file changing the Hamiltonian to do the
 * new model
 * <li> Think which operators you need to construct the Hamiltonian and
 * do the proper updates for them
 * <li> Plot the energies as a function to the magnetic field 
 * (you will to enter a new parameter corresponding to \f$h=\Gamma/J\f$)
 * <li> Plot the entanglement entropy as in the previous exercise. Do you
 * see something when at the QPT?
 * </ul>
 */
