// -*- C++ -*-
/***************************************************************************
 * blitz/rand-normal.h    Random Gaussian (Normal) generator
 *
 * $Id: rand-normal.h,v 1.5 2005/10/19 21:17:00 julianc Exp $
 *
 * Copyright (C) 1997-2001 Todd Veldhuizen <tveldhui@oonumerics.org>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * Suggestions:          blitz-dev@oonumerics.org
 * Bugs:                 blitz-bugs@oonumerics.org
 *
 * For more information, please see the Blitz++ Home Page:
 *    http://oonumerics.org/blitz/
 *
 ***************************************************************************
 *
 * This generator transforms a (0,1] uniform distribution into
 * a Normal distribution.  Let u,v be (0,1] random variables. Then
 *
 *    x = sqrt(-2 ln v) cos(pi (2u-1))
 *
 * is N(0,1) distributed.
 *
 * Reference: Athanasios Papoulis, "Probability, random variables,
 * and stochastic processes," McGraw-Hill : Toronto, 1991.
 *
 ***************************************************************************/

#ifndef BZ_RAND_NORMAL_H
#define BZ_RAND_NORMAL_H

#include <blitz/random.h>
#include <blitz/rand-uniform.h>

#if defined(BZ_HAVE_STD)
#include <cmath>
#else
#include <math.h>
#endif

#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

BZ_NAMESPACE(blitz)

template<typename P_uniform BZ_TEMPLATE_DEFAULT(Uniform)>
class Normal {

public:
    typedef double T_numtype;

    Normal(double mean = 0.0, double variance = 1.0, double = 0.0)
        : mean_(mean), sigma_(BZ_MATHFN_SCOPE(sqrt)(variance))
    { 
    }

    void randomize() 
    { 
        uniform_.randomize();
    }
  
    double random()
    { 
        double u, v;

        do {
            u = uniform_.random();
            v = uniform_.random();    
        } while (v == 0);

        return mean_ + sigma_ * BZ_MATHFN_SCOPE(sqrt)(-2*BZ_MATHFN_SCOPE(log)(v)) * BZ_MATHFN_SCOPE(cos)(M_PI * (2*u - 1));
    } 

private:
    double mean_, sigma_;
    P_uniform uniform_;
};

BZ_NAMESPACE_END

#endif // BZ_RAND_NORMAL_H

