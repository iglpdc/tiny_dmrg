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

blitz::Array<double,2> calculateReducedDensityMatrix(blitz::Array<double,2> psi);
void DMlargeEigen(blitz::Array<double,2>& Hm, blitz::Array<double,2>& Od, 
	const int nn, const int mm);
#endif //HHTRIDI_H 
