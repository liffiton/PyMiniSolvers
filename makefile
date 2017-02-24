CXX ?= g++

r: libminisat.so libminicard.so
d: libminisat.so libminicard.so

r: CFLAGS=-fpic -D NDEBUG -O3 -Wall -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -Wno-parentheses -Wextra
d: CFLAGS=-fpic -D DEBUG -O0 -ggdb -Wall -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -Wno-parentheses -Wextra

OS=$(shell uname -s)
ifeq ($(OS), Darwin)
	SHARED=-dynamiclib
else
	SHARED=-shared
endif

SATINC=minisat/
CARDINC=minicard/

libminisat.so: minisat.o satSolver.o satSystem.o
	$(CXX) $(SHARED) $(CFLAGS) -o $@ $^

minisat.o: minisat.cpp
	$(CXX) -c $(CFLAGS) -I $(SATINC) -o $@ $^

satSolver.o: minisat/minisat/core/Solver.cc
	$(CXX) -c $(CFLAGS) -I $(SATINC) -o $@ $^

satSystem.o: minisat/minisat/utils/System.cc
	$(CXX) -c $(CFLAGS) -I $(SATINC) -o $@ $^
    
libminicard.so: minicard.o cardSolver.o cardSystem.o
	$(CXX) $(SHARED) $(CFLAGS) -o $@ $^

minicard.o: minicard.cpp
	$(CXX) -c $(CFLAGS) -I $(CARDINC) -o $@ $^

cardSolver.o: minicard/minicard/Solver.cc
	$(CXX) -c $(CFLAGS) -I $(CARDINC) -o $@ $^

cardSystem.o: minicard/utils/System.cc
	$(CXX) -c $(CFLAGS) -I $(CARDINC) -o $@ $^

clean:
	rm -f *.so *.o

# check for existence of python versions to control tests
PYTHON2 := $(shell command -v python 2> /dev/null)
PYTHON3 := $(shell command -v python3 2> /dev/null)

test:
ifdef PYTHON2
	python test_minisolvers.py
	python -m doctest -v minisolvers.py
endif
ifdef PYTHON3
	python3 test_minisolvers.py
	python3 -m doctest -v minisolvers.py
endif
	@echo
	@echo "[32mAll tests passed.[m"

.PHONY: clean test
