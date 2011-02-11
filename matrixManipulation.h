#ifndef MATRIX_MANIPULATION_H
#define MATRIX_MANIPULATION_H
 
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

/**
 * @brief A function to reduce a 4 tensor index to a matrix
 *
 * Reduces a m,2,m,2 matrix to a 2m*2m matrix 
 * over the direct product basis of 2 spins
 */
Array<double,2> reduceM2M2(const Array<double,4>& T2, const int m)
{
  int r,c;
  Array<double,2> M4(2*m,2*m);

  r=0;
  for (int a1=0; a1<m; a1++)
    for (int a2=0; a2<2; a2++){
      c=0;
      for (int a3=0; a3<m; a3++)
	for (int a4=0; a4<2; a4++){
	  M4(r,c) = T2(a1,a2,a3,a4);
	  c++;
	}
      r++;    
    }
  return M4;
}
#endif // MATRIX_MANIPULATION_H
