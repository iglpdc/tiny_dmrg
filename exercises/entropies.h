/**
 * @file entropies.h
 * @brief A file with a couple of function to calculate the entanglement
 * entropies 
 * @author Roger Melko 
 * @author Ivan Gonzalez
 * @date $Date$
 *
 * $Revision$ 
 *
 * This is an implementation to give out as the solution of the second
 * exercise in the DMRGloo school.
 *
 * All the functions in this file are inlined. So you can just include
 * them in the original densityMatrix.cpp file. Note that the functions ot
 * calculate the entropies take the reduced denisty matrix eigenvalues,
 * which in the current implementation are only avaliable inside the
 * truncateReducedDM() function.
 *
 * @see truncatedReducedDM()
 */
#ifndef ENTROPIES_H
#define ENTROPIES_H

#include"blitz/array.h"

/**
 * @brief A function to calculate the entanglement entropy 
 *
 * @param density_matrix_eigenvalues the eigenvalues of the reduced
 * density matrix
 *
 * The Von Neumann entanglement entropy is defined as:
 * 
 * \f$S_{vN}(A)\equiv-tr(\rho_{A}\ln\rho_{A})\f$
 *
 * where \f$\rho_{A}\f$ is the reduce density matrix of part \f$A\f$. When
 * \f$\rho_{A}\f$ is diagonalized (\f$\lambda_{i}\f$ being its
 * eigenvalues),
 * you have:
 *
 * \f$S_{vN}(A)=-\sum_{i}\lambda_{i}\ln\lambda_{i}\f$
 *
 * This function assumes that the eigenvalues are OK, i.e.
 * \f$\lambda_{i}\in[0,1],\, \forall i, \quad \sum_{i}\lambda_{i}=1\f$
 */
inline double calculateEntanglementEntropy(const blitz::Array<double,1>& density_matrix_eigenvalues)
{
    double result=0.0;
    for (size_t i=0; i<density_matrix_eigenvalues.size(); ++i)
	result-=density_matrix_eigenvalues(i)*log(density_matrix_eigenvalues(i));
    return result;
    // you can use do it also using the blitz++ way
    // return -1.0*sum(density_matrix_eigenvalues*log(density_matrix_eigenvalues));
}
/**
 * @brief A function to calculate the Renyi entropies 
 *
 * @param density_matrix_eigenvalues the eigenvalues of the reduced
 * density matrix
 *
 * The Renyi entropy of order N is defined as:
 * 
 * \f$S_{n}(A)\equiv \frac{1}{1-n}\ln [tr(\rho^{n}_{A})]\f$
 *
 * where \f$\rho_{A}\f$ is the reduce density matrix of part \f$A\f$. When
 * \f$\rho_{A}\f$ is diagonalized (\f$\lambda_{i}\f$ being its
 * eigenvalues),
 * you have:
 *
 * \f$S_{n}(A)=\frac{1}{1-n}\ln (\sum_{i}\lambda^{n}_{i}})\f$
 *
 * This function assumes that the eigenvalues are OK, i.e.
 * \f$\lambda_{i}\in[0,1],\, \forall i, \quad \sum_{i}\lambda_{i}=1\f$

 */
template<int N>
inline double calculateRenyiEntropy(const blitz::Array<double,1>& density_matrix_eigenvalues)
{
    double result=0.0;
    for (size_t i=0; i<density_matrix_eigenvalues.size(); ++i)
	result=pow(density_matrix_eigenvalues(i),N);
    return log(result)/(1-N);
    // you can use do it also using the blitz++ way
    // return log(sum(pow(density_matrix_eigenvalues,N)))/(1-N);
}
/**
 * @brief A specialization of function to calculate the Renyi entropy for N=1 
 *
 * @param density_matrix_eigenvalues the eigenvalues of the reduced
 * density matrix
 * 
 * \f$\lim_{n\to 1}S_{n}(A)=S_{vN}(A)
 *
 * This specialization avoids the problem of dividing by zero.
 */
template<>
inline double calculateRenyiEntropy<1>(const blitz::Array<double,1>& density_matrix_eigenvalues)
{
    return calculateEntanglementEntropy(density_matrix_eigenvalues);
}
#endif // ENTROPIES_H
