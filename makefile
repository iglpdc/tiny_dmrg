OBJS = Heis_FSA.o HHtridi8.o lanczosDMRG.o

a.out: $(OBJS)
	g++ -O3 $(OBJS)
lanczosDMRG.o: lanczosDMRG.cpp lanczosDMRG.h
	g++ -c -O3 lanczosDMRG.cpp
HHtridi8.o: HHtridi8.cpp heis_dmrg.h
	g++ -c -O3 HHtridi8.cpp
Heis_FSA.o: Heis_FSA.cpp heis_dmrg.h
	g++ -c -O3 Heis_FSA.cpp

clean :
	rm *.o
