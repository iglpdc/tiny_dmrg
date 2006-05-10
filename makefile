a.out: HeisFSA_2.o HHtridi8.o LanczosDMRG_4.o 
	g++ HeisFSA_2.o HHtridi8.o LanczosDMRG_4.o
LanczosDMRG_.o: LanczosDMRG_4.cpp
	g++ -c -O3 LanczosDMRG_4.cpp
HHtridi8.o: HHtridi8.cpp
	g++ -c -O3 HHtridi8.cpp
HeisFSA_2.o: HeisFSA_2.cpp
	g++ -c -O3 HeisFSA_2.cpp
