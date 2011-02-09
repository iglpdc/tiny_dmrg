#ifndef LANCZOS_DMRG_H
#define LANCZOS_DMRG_H
 
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int LanczosED(Array<double,2>&, Array<double,1>&, double *, const int);
void Normalize(Array<double,1>& V);
int tqli2(Array<double,1>& d, Array<double,1>& e, int n, Array<double,2>& z, 
	const int Evects);
#endif
