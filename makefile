OBJS = Heis_FSA.o HHtridi8.o LanczosDMRG_4.o

a.out: $(OBJS)
	g++ $(OBJS)
LanczosDMRG_4.o: LanczosDMRG_4.cpp heis_dmrg.h
	g++ -c -O LanczosDMRG_4.cpp
HHtridi8.o: HHtridi8.cpp heis_dmrg.h
	g++ -c -O HHtridi8.cpp
Heis_ FSA.o: Heis_FSA.cpp heis_dmrg.h
	g++ -c -O Heis_FSA.cpp

clean :
	rm *.o
