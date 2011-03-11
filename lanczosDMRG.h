/**
 * @file lanczosDMRG.h
 *
 * @brief headers for Lanczos ED
 *
 * @author Roger Melko 
 * @author Ivan Gonzalez
 * $Date$
 *
 * $Revision$ 
 */
#ifndef LANCZOS_DMRG_H
#define LANCZOS_DMRG_H

#include"blitz/array.h"
 
double calculateGroundState(blitz::Array<double,4>&, blitz::Array<double,2>&);
int diagonalizeWithLanczos(blitz::Array<double,2>&, blitz::Array<double,1>&, double *);
#endif // LANCZOS_DMRG_H
