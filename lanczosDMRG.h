#ifndef LANCZOS_DMRG_H
#define LANCZOS_DMRG_H
 
double calculateGroundState(blitz::Array<double,4>&, blitz::Array<double,2>&, 
	const int);
int LanczosED(blitz::Array<double,2>&, blitz::Array<double,1>&, double *, 
	const int);
#endif // LANCZOS_DMRG_H
