/** 
 * @file HHtridi8.h
 *
 * @brief Interface for the routines related to the calculation of the  
 * reduced density matrix
 *
 * Oct. 26, Roger Melko ORNL
 */
#ifndef HHTRIDI_H
#define HHTRIDI_H  

blitz::Array<double,2> transformOperator(const blitz::Array<double,2>& op, 
	const blitz::Array<double,2>& transposed_transformation_matrix,
	const blitz::Array<double,2>& transformation_matrix);

blitz::Array<double,2> calculateReducedDensityMatrix(blitz::Array<double,2> psi);

blitz::Array<double,2> truncateReducedDM(blitz::Array<double,2>& Hm, 
	const int mm);
#endif //HHTRIDI_H 
