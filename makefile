
r: libminisat.so libminicard.so
d: libminisat.so libminicard.so

r: CFLAGS=-fpic -D NDEBUG -O3 -Wall -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS
d: CFLAGS=-fpic -D DEBUG -O0 -ggdb -Wall -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS

SATINC=minisat/
CARDINC=minicard/

libminisat.so: minisat.o satSolver.o satSystem.o
	g++ -shared $(CFLAGS) -o $@ $^

minisat.o: minisat.cpp
	g++ -c $(CFLAGS) -I $(SATINC) -o $@ $^

satSolver.o: minisat/minisat/core/Solver.cc
	g++ -c $(CFLAGS) -I $(SATINC) -o $@ $^

satSystem.o: minisat/minisat/utils/System.cc
	g++ -c $(CFLAGS) -I $(SATINC) -o $@ $^
    
libminicard.so: minicard.o cardSolver.o cardSystem.o
	g++ -shared $(CFLAGS) -o $@ $^

minicard.o: minicard.cpp
	g++ -c $(CFLAGS) -I $(CARDINC) -o $@ $^

cardSolver.o: minicard/minicard/Solver.cc
	g++ -c $(CFLAGS) -I $(CARDINC) -o $@ $^

cardSystem.o: minicard/utils/System.cc
	g++ -c $(CFLAGS) -I $(CARDINC) -o $@ $^

clean:
	rm -f *.so *.o

.PHONY: clean
