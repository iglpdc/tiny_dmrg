/** 
 * @file HHtridi8.cpp
 *
 * @brief Implementation for the routines related to the calculation of the  
 * reduced density matrix
 *
 * Oct. 26, Roger Melko ORNL
 */

#include "blitz/array.h"
#include "exceptions.h"
#include "tred3.h"
#include "tqli2.h"
#include "HHtridi8.h"

/**
 * @brief A function to calculate the reduced density matrix 
 *
 * @param Psi the wavefunction with you want to calculate the density matrix 
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
 * @param truncated_density_matrix the truncated reduced density matrix
 * @param nn total number of density matrix eigenvalues
 * @param mm number of density matrix eigenvalues to keep
 *
 * Takes the reduced density matrix (DM) which is a real and symmetric matrix. 
 * Performs a Householder reduction to tridiagonal form, then diagonalizes
 * exactly this matrix. Takes input as Blitz++ arrays 
 *
 * checks if the DM is square (but not if it's symmetric)
 * checks that mm=<nn
 * assures that the sum of the DM eigenvalues is 1.0
 *
 */
void DMlargeEigen(blitz::Array<double,2>& density_matrix, 
	blitz::Array<double,2>& truncated_density_matrix, 
	const int nn, const int mm)
{
    if (density_matrix.cols()!=density_matrix.rows())
	throw dmrg::Exception("reduced DM is not square");
    
    if (mm>nn)
	throw dmrg::Exception("Cannot keep more states than size of DM");

    blitz::Array<double,1> e(nn);
    blitz::Array<double,1> density_matrix_eigenvalues(nn); 

    // Householder reduction: reduces symmetric matrix to a tridiagonal
    // form
    tred3(density_matrix, density_matrix_eigenvalues, e, nn);

    // diagonalizes a triangular matrix
    int rtn = tqli2(density_matrix_eigenvalues, e, nn, density_matrix, 1);

    // now density_matrix(j,i) is the eigenvector corresponding to d[i]

    // check that the sum of the eigenvalues is 1.0 and stop if not
    double sum_of_density_matrix_eigenvalues=sum(density_matrix_eigenvalues);

    if (fabs(1.0-sum_of_density_matrix_eigenvalues) > 0.00001)
	std::cout<<"sum_of_density_matrix_eigenvalues error:"\
	    <<sum_of_density_matrix_eigenvalues<<std::endl;	

    blitz::Array<int,1> inx(nn);

    for (int j=0; j<nn; j++) inx(j)=j;

    //STRIAGHT INSERTION SORT O(N^2)    
    for (int j=1; j<nn; j++)
    {         
	double a = density_matrix_eigenvalues(j);
	int b = inx(j);
	int i=j-1;
	while (i>=0 && density_matrix_eigenvalues(i) >a)
	{
	    density_matrix_eigenvalues(i+1)=density_matrix_eigenvalues(i);
	    inx(i+1)=inx(i);
	    i--;
	}
	density_matrix_eigenvalues(i+1)=a;
	inx(i+1)=b;
    }
    //std::cout<<"density_matrix_eigenvalues\n"<<density_matrix_eigenvalues;
    //std::cout<<"\n inx\n"<<inx;

    // calculate the truncated error
    double truncation_error = 0;
    for (int kk=nn-1; kk<(nn-1-mm); kk--)
    {
	truncation_error += density_matrix_eigenvalues(kk);
    }
    //std::cout<<"truncation_error "<<" "<<truncation_error<<std::endl;

    // define the truncation matrix formed by the eigenvector corresponding 
    // to the largest eigevalues
    for (int  kk=0; kk<mm; kk++)
    {
	for (int i=0; i<nn; i++)
	{
	    truncated_density_matrix(kk,i) = density_matrix(i,inx(nn-1-kk));   
	}
    } 
}
// end HHtridi8.cpp
