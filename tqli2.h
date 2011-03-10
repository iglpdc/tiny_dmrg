/**
 * @file tqli2.h
 *
 * @brief Interface for the tqli2 function using blitz arrays
 *
 * @author Roger Melko 
 * @author Ivan Gonzalez
 * @date $Date$
 *
 * $Revision$ 
 */
#ifndef TQLI2_H
#define TQLI2_H

#include"blitz/array.h"

int tqli2(blitz::Array<double,1>&, blitz::Array<double,1>&, int, blitz::Array<double,2>& , const int);
#endif // TQLI2_H
