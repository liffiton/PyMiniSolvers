
CFLAGS=-fpic -D NDEBUG -O3 -Wall -I minisat/ -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS

libminisat.so: minisat.o minisat/minisat/core/Solver.o minisat/minisat/utils/System.o
	g++ -shared $(CFLAGS) -o $@ $^

minisat.o: minisat.cpp
	g++ -c $(CFLAGS) -o $@ $^

minisat/minisat/%.o: minisat/minisat/%.cc
	g++ -c $(CFLAGS) -o $@ $^

clean:
	rm -f libminisat.so minisat.o minisat/minisat/core/Solver.o minisat/minisat/utils/System.o

.PHONY: clean
