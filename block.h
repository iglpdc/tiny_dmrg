#ifndef BLOCK_H 
#define BLOCK_H

#include <blitz/array.h>
#include <fstream>

BZ_USING_NAMESPACE(blitz)

//block class
class BLOCK {
	public:
		// number of sites in the block
		int size;    
		///A' plus right spin Hamiltonian
		Array<double,2> HAB;   

		BLOCK();
		void ISAwrite(const int sites);
		void FSAread(const int sites,const int iter);
		void FSAwrite(const int sites,const int iter);

	private:
	    //filename for storing block on disk
		char fname[7];
		void Read();
		void Write();
};
BLOCK::BLOCK(){
///constructor initialize filename
  fname[0] = '.'; fname[1] = 48;  //ASCII for 0
  fname[4] = '.'; fname[5] = 'r'; 
  fname[6] = '\0';
}

void BLOCK::ISAwrite(const int sites){
	/// rename/write wrapper for the ISA
	fname[3] = 48 + (sites)%10;          //some ASCII crap
	fname[2] = 48 + sites/10;
	fname[5] = 'l';
	Write(); //write left block
	fname[5] = 'r'; 
	Write(); //write right block
}
void BLOCK::FSAread(const int sites,const int iter){
	fname[3] = 48 + (sites)%10;          //some ASCII crap
	fname[2] = 48 + sites/10;
	if (iter%2 == 0) fname[5] = 'r';  else fname[5]= 'l';
	Read();
}//FSAread
void BLOCK::FSAwrite(const int sites,const int iter){
      fname[3] = 48 + (sites)%10;          //some ASCII crap
      fname[2] = 48 + sites/10;
      if (iter%2 == 0) fname[5] = 'l';  else fname[5]= 'r';
      Write();
}//FSAwrite

void BLOCK::Write()
{
  ofstream fout;  
  fout.open(fname,ios::out);
  //  fout << blk->size <<endl;
  fout <<setprecision(12)<<HAB ;
  fout.close();

} //Write
void BLOCK::Read()
{
  ifstream fin;  
  fin.open(fname,ios::in);
  //  fin >> blk->size; 
  fin >> HAB ;
  fin.close();
//  cout<<blk->HAp;
}//Read

#endif
