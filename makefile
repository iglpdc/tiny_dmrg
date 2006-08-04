OBJS = HeisFSA_3.o HHtridi8.o LanczosDMRG_4.o

a.out: $(OBJS)
	g++ $(OBJS)
LanczosDMRG_4.o: LanczosDMRG_4.cpp heis_dmrg.h
	g++ -c -O LanczosDMRG_4.cpp
HHtridi8.o: HHtridi8.cpp heis_dmrg.h
	g++ -c -O HHtridi8.cpp
HeisFSA_3.o: HeisFSA_3.cpp heis_dmrg.h
	g++ -c -O HeisFSA_3.cpp

clean :
	rm *.o
