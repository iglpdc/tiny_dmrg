a.out: HeisInf_6.o HHtridi8.o LanczosDMRG_4.o 
	g++ HeisInf_6.o HHtridi8.o LanczosDMRG_4.o
LanczosDMRG_.o: LanczosDMRG_4.cpp
	g++ -c -O3 LanczosDMRG_4.cpp
HHtridi8.o: HHtridi8.cpp
	g++ -c -O3 HHtridi8.cpp
HeisInf_6.o: HeisInf_6.cpp
	g++ -c -O3 HeisInf_6.cpp
