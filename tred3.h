//Oct. 26, Roger Melko ORNL
//
//Exact diagonalization program for a real, symmetric matrix
//Householder reduction to tridiagonal form, then ED
//
//Takes input as Blitz++ Arrays<>
//
//set the below value to "1" if you want the eigenvectors, and to "0" if not
#ifndef TRED3_H
#define TRED3_H

void tred3(blitz::Array<double,2>& a, blitz::Array<double,1>& d, blitz::Array<double,1>& e, int n);
#endif // TRED3_H
