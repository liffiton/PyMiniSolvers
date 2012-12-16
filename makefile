
CFLAGS=-fpic -D NDEBUG -O3 -Wall -I minisat/ -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS

libminisat.so: minisat.o Solver.o System.o
	g++ -shared $(CFLAGS) -o $@ $^

minisat.o: minisat.cpp
	g++ -c $(CFLAGS) -o $@ $^

Solver.o: minisat/minisat/core/Solver.cc
	g++ -c $(CFLAGS) -o $@ $^

System.o: minisat/minisat/utils/System.cc
	g++ -c $(CFLAGS) -o $@ $^
    
clean:
	rm -f libminisat.so *.o

.PHONY: clean
