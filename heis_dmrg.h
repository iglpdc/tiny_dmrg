#ifndef HEIS_DMRG_H
#define HEIS_DMRG_H
 
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

//block class
	class BLOCK {
		public:
		    // number of sites in the block
			int size;    

			///A' plus right spin Hamiltonian
			Array<double,2> HAB;   

		private:
	};


//template function prototypes
int tqli2(Array<double,1>&, Array<double,1>&, int, Array<double,2>&, const int);
void tred3(Array<double,2>& , Array<double,1>& , Array<double,1>& , int );

template<typename T> Array<T,2> reduceM2M2(const Array<T,4>&, const int);
void EigenValuesLAN(Array<double,4>&, Array<double,2>&, const int, double *);
void DMlargeEigen(Array<double,2>&, Array<double,2>&, const int, const int);
void BlockWrite(BLOCK *, const int, char []);
void BlockRead(BLOCK *, const int, char []);


#endif
