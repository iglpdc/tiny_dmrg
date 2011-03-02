OPT:=3

default: incremental
exec := a.out
rev:= $(shell svnversion -n)

ifdef debug 
CXXFLAGS +=-O$(OPT) -pg -Wall -I.
else
CXXFLAGS +=-O$(OPT) -I.
endif

ifdef threads 
CXXFLAGS+=-pthreads
endif

OBJS = tqli2.o tred3.o Heis_FSA.o densityMatrix.o lanczosDMRG.o

$(exec): $(OBJS)
	g++ $(CXXFLAGS) $(OBJS)
tqli2.o: tqli2.cpp tqli2.h
	g++ -c $(CXXFLAGS) tqli2.cpp
tred3.o: tred3.cpp tred3.h
	g++ -c $(CXXFLAGS) tred3.cpp
lanczosDMRG.o: lanczosDMRG.cpp lanczosDMRG.h lanczosDMRG_helpers.h tqli2.o
	g++ -c $(CXXFLAGS) lanczosDMRG.cpp
densityMatrix.o: densityMatrix.cpp tred3.o tqli2.o
	g++ -c $(CXXFLAGS) densityMatrix.cpp
Heis_FSA.o: Heis_FSA.cpp matrixManipulation.h
	g++ -c $(CXXFLAGS) Heis_FSA.cpp

.PHONY: clean incremental all doc tarball

all: clean incremental doc

clean:
	rm *.o

doc:
	doxygen ./bin/doxy_sse.cfg

tarball:
	./bin/make_tarball

incremental: ${exec}
