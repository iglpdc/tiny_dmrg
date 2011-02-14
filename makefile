OPT:=3

default: incremental
exec := a.out

ifdef debug 
CXXFLAGS+=-O$(OPT) -gp -Wall
else
CXXFLAGS+=-O$(OPT) -Wall
endif

ifdef threads 
CXXFLAGS+=-pthreads
endif

OBJS = tqli2.o tred3.o Heis_FSA.o HHtridi8.o lanczosDMRG.o

$(exec): $(OBJS)
	g++ -O3 $(OBJS)
tqli2.o: tqli2.cpp tqli2.h
	g++ -c -O3 tqli2.cpp
tred3.o: tred3.cpp tred3.h
	g++ -c -O3 tred3.cpp
lanczosDMRG.o: lanczosDMRG.cpp lanczosDMRG.h lanczosDMRG_helpers.h tqli2.o
	g++ -c -O3 lanczosDMRG.cpp
HHtridi8.o: HHtridi8.cpp tred3.o tqli2.o
	g++ -c -O3 HHtridi8.cpp
Heis_FSA.o: Heis_FSA.cpp matrixManipulation.h
	g++ -c -O3 Heis_FSA.cpp

.PHONY: clean tarball incremental all

all: clean incremental

clean :
	rm *.o

incremental: ${exec}
