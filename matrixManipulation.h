#ifndef MATRIX_MANIPULATION_H
#define MATRIX_MANIPULATION_H
 
#include "blitz/array.h"

/**
 * @brief A function to reduce a 4-index tensor to a matrix
 *
 * Reduces a m,2,m,2 matrix to a 2m*2m matrix 
 * over the direct product basis of 2 spins
 */
inline blitz::Array<double,2> reduceM2M2(const blitz::Array<double,4>& tensor, const int firstDim, const int secondDim)
{
    const int matrixDim=firstDim*secondDim;
      
    blitz::Array<double,2> result(matrixDim, matrixDim);

    int c;
    int r=0;
    for (int a1=0; a1<firstDim; a1++)
    {
      for (int a2=0; a2<secondDim; a2++)
      {
	  c=0;
	  for (int a3=0; a3<firstDim; a3++)
	  {
	      for (int a4=0; a4<secondDim; a4++)
	      {
		  result(r,c) = tensor(a1,a2,a3,a4);
		  c++;
	      }
	  }
	  r++;    
      }
    }
    return result;
}
#endif // MATRIX_MANIPULATION_H
