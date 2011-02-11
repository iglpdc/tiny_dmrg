#ifndef LANCZOS_DMRG_HELPERS_H
#define LANCZOS_DMRG_HELPERS_H
 
#include <blitz/array.h>
#include <cmath>

BZ_USING_NAMESPACE(blitz)

/**
 * @brief A function to do the dot prodcut of two wavefunctions
 *
 * @param V1 the wavefunction to multiply
 * @param V2 the wavefunction to multiply
 */
double dotProduct(const Array<double,1>& V1, const Array<double,1>& V2)
{
    return sum(V1*V2);
}
/**
 * @brief A function to get the norm of a wavefunction
 *
 * @param V the wavefunction to normalize
 *
 * Takes a wavefunction and return its norm, i.e. the square 
 * root of dot product with itself
 */
double calculateNorm(const Array<double,1>& V) 
{
  double norm = dotProduct(V,V);           
  return sqrt(norm);
}

/**
 * @brief A function to randomize a wavefunction
 *
 * @param V the wavefunction to randomize
 *
 * Takes a wavefunction fills all its components with a random real number
 */
void randomize(Array<double,1>& V) 
{
  for (int i=0; i<V.size(); i++)
  {
      V(i) = rand()%10*0.1;  //random starting vec
      if ( (rand()%2) == 0) V(i) *= -1.0000001;
  }
}

/**
 * @brief A function to normalize a wavefunction
 *
 * @param V the wavefunction to normalize
 *
 * Takes a wavefunction and normalizes it
 */
void Normalize(Array<double,1>& V) 
{
  double norm = calculateNorm(V);

  for (int i=0; i<V.size(); i++)
      V(i) /= norm;
}
#endif // LANCZOS_DMRG_HELPERS_H
