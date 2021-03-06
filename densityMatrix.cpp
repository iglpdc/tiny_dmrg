/** 
 * @file densityMatrix.cpp
 *
 * @brief Implementation for the routines related to the calculation of the  
 * reduced density matrix
 *
 * @author Roger Melko 
 * @author Ivan Gonzalez
 * @date $Date$
 *
 * $Revision$ 
 */
#include "blitz/array.h"
#include "exceptions.h"
#include "tred3.h"
#include "tqli2.h"
#include "densityMatrix.h"

/**
 * @brief A function to transform an operator to the new (truncated) basis
 *
 * @param op matrix with the operator in the old basis
 * @param transformation_matrix a matrix transforming the old basis to the new 
 * (truncated) one 
 * @param transposed_transformation_matrix the transpose of the
 * transformation matrix
 *
 * @return a matrix with the transformed operator
 *
 * The transformation matrix should be the (truncated) reduced density
 * matrix that you obtain as a result of applying the
 * truncateDensityMatrix function
 */
blitz::Array<double,2> transformOperator(const blitz::Array<double,2>& op, 
	const blitz::Array<double,2>& transposed_transformation_matrix,
	const blitz::Array<double,2>& transformation_matrix)
{
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::thirdIndex k;

    blitz::Array<double,2> tmp(op.rows(),
	    transposed_transformation_matrix.cols());
    tmp=sum(op(i,k)*transposed_transformation_matrix(k,j),k);

    blitz::Array<double,2> result(transformation_matrix.rows(), 
	    tmp.cols());
    result=sum(transformation_matrix(i,k)*tmp(k,j),k);

    return result;
}

/**
 * @brief A function to calculate the reduced density matrix 
 *
 * @param psi the wavefunction with you want to calculate the density matrix 
 *
 * @return a matrix with the reduced density matrix
 *
 * The wavefunction has to be written as a matrix.
 *
 */
blitz::Array<double,2> calculateReducedDensityMatrix(blitz::Array<double,2> psi)
{
    int rows_psi=psi.rows();
    int cols_psi=psi.cols();

    blitz::Array<double,2> result(rows_psi, rows_psi);
    result=0.0;

    for (int i=0; i<rows_psi; i++)
      for (int j=0; j<rows_psi; j++)
        for (int k=0; k<cols_psi; k++)
 	  result(i,j) += psi(i,k)*psi(j,k); 

    return result;
}

/**
 * @brief A function calculate the density matrix eigenvalues
 *
 * @param density_matrix the reduced density matrix
 * @param m number of density matrix eigenvalues to keep
 *
 * @return truncated_density_matrix the truncated reduced density matrix
 *
 * Takes the reduced density matrix (DM) which is a real and symmetric 
 * matrix. Performs a Householder reduction to tridiagonal form, then 
 * diagonalizes exactly this matrix. Takes input as Blitz++ arrays 
 *
 * Checks if the DM is square (but not if it's symmetric), and that
 * mm=<nn. Assures that the sum of the DM eigenvalues is 1.0
 *
 */
blitz::Array<double,2> truncateReducedDM(blitz::Array<double,2>& 
	density_matrix, int m)
{
    if (density_matrix.cols()!=density_matrix.rows())
	throw dmrg::Exception("reduced DM is not square");
    
    const int n=density_matrix.rows();

    if (m>n)
	throw dmrg::Exception("Cannot keep more states than size of DM");

    blitz::Array<double,1> density_matrix_eigenvalues(n); 

    diagonalizeDensityMatrix(density_matrix, density_matrix_eigenvalues);

    // check that the sum of the eigenvalues is close to 1.0
    if (fabs(1.0-sum(density_matrix_eigenvalues)) > 0.00001)
	throw dmrg::Exception("sum_of_density_matrix_eigenvalues is not one");

    // get the indexes of the largest eigenvalues
    blitz::Array<int,1> indexes=
	orderDensityMatrixEigenvalues(density_matrix_eigenvalues);

    // calculate and print the truncation error
    //std::cout<<std::setprecision(16)<<"truncation_error ";
    //std::cout<<calculateTruncationError(density_matrix_eigenvalues, indexes, mm)<<'\n';

    // define the truncation matrix formed by the eigenvector corresponding 
    // to the largest eigevalues

    blitz::Array<double,2> truncated_density_matrix(m, n); 

    for (int  kk=0; kk<m; kk++)
    {
	for (int i=0; i<n; i++)
	{
	    truncated_density_matrix(kk,i) = density_matrix(i,indexes(n-1-kk));   
	}
    }
    return truncated_density_matrix; 
}
/** 
 * @brief A function to order the reduced density matrix eigenvalues
 *
 * Given the (reduced) density matrix eigenvalues in an array this function
 * returns an array with the permutaion of the indexes corresponding to
 * ordering the eigenvalues in decreasing order.
 *
 * @param density_matrix_eigenvalues the (reduced) density matrix
 * eigenvalues
 *
 * @result an array with the permutation of the indexes corresponding to
 * ordering the eigenvalues in decreasing order. 
 *
 * E.g. if density_matrix_eigenvalues={ 0.15, 0.8, 0.05 }, this function returns
 * {1,0,2}
 *
 */
blitz::Array<int,1> orderDensityMatrixEigenvalues(
	blitz::Array<double,1>& density_matrix_eigenvalues) 
{ 
    int nn=density_matrix_eigenvalues.size(); 
    
    if (nn<1)
	throw dmrg::Exception("orderDensityMatrixEigenvalues: no \
		eigenvalues passed");

    blitz::Array<int,1> result(nn);
    for (int j=0; j<nn; j++) result(j)=j;

    //STRIAGHT INSERTION SORT O(N^2)    
    for (int j=1; j<nn; j++)
    {         
	double a = density_matrix_eigenvalues(j);
	int b = result(j);
	int i=j-1;
	while (i>=0 && density_matrix_eigenvalues(i) >a)
	{
	    //density_matrix_eigenvalues(i+1)=density_matrix_eigenvalues(i);
	    result(i+1)=result(i);
	    i--;
	}
	//density_matrix_eigenvalues(i+1)=a;
	result(i+1)=b;
    }
    return result;
}
/**
 * @brief A function to diagonalize the reduced density matrix
 *
 * @param density_matrix a matrix with the reduced density matrix
 * @param density_matrix_eigenvalues an array with the reduced density
 * matrix eigenvalues
 * 
 * On entrace to the function, density_matrix is the reduced density
 * matrix, and density_matrix_eigenvalues can be garbage.
 * On return, density_matrix(j,i) is the eigenvector corresponding to 
 * density_matrix_eigenvalues(i). To diagonalize the density matrix, which
 * is symmetric, you first reduce it to a tridiagonal form using
 * Householder reduction, and then diagonalize this tridiagonal matrix. 
 */
void diagonalizeDensityMatrix(blitz::Array<double,2>& 
	density_matrix, blitz::Array<double,1>& density_matrix_eigenvalues)
{
    const int n=density_matrix.rows();
    // temporary array
    blitz::Array<double,1> e(n);

    // reduce symmetric matrix to a tridiagonal form (Householder reduct.)
    tred3(density_matrix, density_matrix_eigenvalues, e, n);

    // diagonalizes a tridiagonal matrix
    int rtn = tqli2(density_matrix_eigenvalues, e, n, density_matrix, 1);
} 
/**
 * @brief A function to calculate the truncation error
 *
 * @param density_matrix_eigenvalues all of them
 * @param indexes the permutation of the indexes of the eigenvalues that
 * orders them in increasing value
 * @param m number of states to keep
 *
 * @return a double with the truncation error
 *
 * The truncation error is given by the sum of the reduced density matrix
 * eigenvalues corresponding to the eigenvectors that are truncated out
 *
 */
double calculateTruncationError(const blitz::Array<double,1>& density_matrix_eigenvalues, 
	const blitz::Array<int,1>& indexes, int m)
{
    double result = 0.0;
    for (size_t i=0; i<density_matrix_eigenvalues.size()-m; ++i)
    {
	result+=density_matrix_eigenvalues(indexes(i));
    }
    return result;
}
// end densityMatrix.cpp
