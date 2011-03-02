/** 
 * @file densityMatrix.h
 *
 * @brief Interface for the routines related to the calculation of the  
 * reduced density matrix
 *
 * Oct. 26, Roger Melko ORNL
 */
#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H  

#include "blitz/array.h"

blitz::Array<double,2> transformOperator(const blitz::Array<double,2>& op, 
	const blitz::Array<double,2>& transposed_transformation_matrix,
	const blitz::Array<double,2>& transformation_matrix);

blitz::Array<double,2> calculateReducedDensityMatrix(blitz::Array<double,2> psi);

blitz::Array<double,2> truncateReducedDM(blitz::Array<double,2>& density_matrix, 
	const int mm);

blitz::Array<int,1> orderDensityMatrixEigenvalues(
	blitz::Array<double,1>& density_matrix_eigenvalues);

double calculateTruncationError(
	const blitz::Array<double,1>& density_matrix_eigenvalues, 
	const blitz::Array<int,1>& indexes, int m);

#endif //DENSITY_MATRIX_H 
