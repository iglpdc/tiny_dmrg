/**
 * @file block.h
 *
 * @brief A file that contains the block class used in the DMRG
 *
 * @author Roger Melko 
 * @author Ivan Gonzalez
 * $Date$
 *
 * $Revision$ 
 */
#ifndef BLOCK_H 
#define BLOCK_H

#include <fstream>
#include "blitz/array.h"

//BZ_USING_NAMESPACE(blitz)

///Block class
class Block {
	public:
		/// number of sites in the block
		int size;    
		/// A' plus right spin Hamiltonian: Blitz++ array
		blitz::Array<double,2> HAB;   

		Block();
		void ISAwrite(const int sites);
		void FSAread(const int sites,const int iter);
		void FSAwrite(const int sites,const int iter);

	private:
	    ///filename for storing block on disk
		char fname[7];

		void Read();
		void Write();
};
Block::Block(){
///constructor: initialize filename for writing blocks to file
  fname[0] = '.'; fname[1] = 48;  //ASCII for 0
  fname[4] = '.'; fname[5] = 'r'; 
  fname[6] = '\0';
}

void Block::ISAwrite(const int sites){
/// rename/write wrapper for the ISA
	fname[3] = 48 + (sites)%10;          //some ASCII 
	fname[2] = 48 + sites/10;
	fname[5] = 'l';
	Write(); //write left block
	fname[5] = 'r'; 
	Write(); //write right block
}

void Block::FSAread(const int sites,const int iter){
/// file read for the finite-system algorithm
	fname[3] = 48 + (sites)%10;          //some ASCII 
	fname[2] = 48 + sites/10;
	if (iter%2 == 0) fname[5] = 'r';  else fname[5]= 'l';
	Read();
}//FSAread

void Block::FSAwrite(const int sites,const int iter){
/// file write for the finite-system algorithm
      fname[3] = 48 + (sites)%10;         //some ASCII 
      fname[2] = 48 + sites/10;
      if (iter%2 == 0) fname[5] = 'l';  else fname[5]= 'r';
      Write();
}//FSAwrite

void Block::Write() {
/// opens the output file and writes the Blitz++ array
  std::ofstream fout;  
  fout.open(fname,std::ios::out);
  fout <<std::setprecision(12)<<HAB ;
  fout.close();
} //Write

void Block::Read() {
/// opens the input file and reads the Blitz++ array
  std::ifstream fin;  
  fin.open(fname,std::ios::in);
  fin >> HAB ;
  fin.close();
}//Read

#endif
