/** 
 * @file matrixManipulation.h
 * @brief Functions for handling matrices
 * 
 * @author Roger Melko 
 * @author Ivan Gonzalez
 * @date $Date$
 *
 * $Revision$ 
 */
#ifndef MATRIX_MANIPULATION_H
#define MATRIX_MANIPULATION_H
 
#include "blitz/array.h"
#include "exceptions.h"

/**
 * @brief A function to create a identity matrix
 */
inline blitz::Array<double,2> createIdentityMatrix(int rows)
{
    if (!(rows>0))
	throw dmrg::Exception("createIdentity: rows -ve");
    blitz::Array<double,2> result(rows,rows);
    result=0.0;
    for (int i=0; i<rows; ++i)
	result(i,i)=1.0;
    return result;
}

/**
 * @brief A function to reduce a 4-index tensor to a matrix
 *
 * Reduces a m,2,m,2 matrix to a 2m*2m matrix 
 * over the direct product basis of 2 spins
 */
inline blitz::Array<double,2> reduceM2M2(const blitz::Array<double,4>& tensor)
{
    const int first_dim=tensor.extent(blitz::firstDim);
    const int second_dim=tensor.extent(blitz::secondDim);

    if (first_dim!=tensor.extent(blitz::thirdDim))
	throw dmrg::Exception("reduceTensorToMatrix: wrong dims");

    if (second_dim!=tensor.extent(blitz::fourthDim))
	throw dmrg::Exception("reduceTensorToMatrix: wrong dims");

    const int matrixDim=first_dim*second_dim;
      
    blitz::Array<double,2> result(matrixDim, matrixDim);

    int c;
    int r=0;
    for (int a1=0; a1<first_dim; a1++)
    {
      for (int a2=0; a2<second_dim; a2++)
      {
	  c=0;
	  for (int a3=0; a3<first_dim; a3++)
	  {
	      for (int a4=0; a4<second_dim; a4++)
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
