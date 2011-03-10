/**
 * @file tred3.h
 *
 * @brief Interface for the tred3 function using blitz arrays
 *
 * @author Roger Melko 
 * @author Ivan Gonzalez
 * @date $Date$
 *
 * $Revision$ 
 */
#ifndef TRED3_H
#define TRED3_H

#include "blitz/array.h"

void tred3(blitz::Array<double,2>& a, blitz::Array<double,1>& d, blitz::Array<double,1>& e, int n);
#endif // TRED3_H
