//Oct. 26, Roger Melko ORNL
//
//Exact diagonalization program for a real, symmetric matrix
//Householder reduction to tridiagonal form, then ED
//
//Takes input as Blitz++ Arrays<>
//
//set the below value to "1" if you want the eigenvectors, and to "0" if not
#ifndef HHTRIDI_H
#define HHTRIDI_H  

void DMlargeEigen(blitz::Array<double,2>& Hm, blitz::Array<double,2>& Od, const int nn, 
		  const int mm);
#endif //HHTRIDI_H 
