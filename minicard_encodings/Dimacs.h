/****************************************************************************************[Dimacs.h]
MiniCARD Copyright (c) 2012, Mark Liffiton, Jordyn Maglalang

MiniCARD is based on MiniSAT, whose original copyright notice is maintained below,
and it is released under the same license.

---

Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef Minisat_Dimacs_h
#define Minisat_Dimacs_h

#include <stdio.h>

#include "utils/ParseUtils.h"
#include "core/SolverTypes.h"

namespace Minisat {

//=================================================================================================
// DIMACS Parser:

// readConstr returns :
//      false if a clause was found
//      true if an atmost was found
//      true if an atleast was found (and modifies the constriant to be an atmost)
template<class B, class Solver>
static int readConstr(B& in, Solver& S, vec<Lit>& lits, int& bound) {
    int     parsed_lit, var;
    lits.clear();
    bound = 0;
    for (;;){
        skipWhitespace(in);
        if (eagerMatch(in,"<=")) goto AtMost;
        if (eagerMatch(in,">=")) goto AtLeast;
        parsed_lit = parseInt(in);
        if (parsed_lit == 0) break;
        var = abs(parsed_lit)-1; 
        lits.push( (parsed_lit > 0) ? mkLit(var) : ~mkLit(var) );
    }
    return false;
    
    AtMost:
    bound = parseInt(in);
    return true;
    
    AtLeast:
    bound = lits.size() - parseInt(in);
    for (int i=0;i<lits.size();i++) {
        lits[i] = ~lits[i];
    }
    return true;
}

template<class B, class Solver>
static void parse_DIMACS_main(B& in, Solver& S) {
    vec<Lit> lits;
    int bound = 0;
    
    int vars    = 0;
    int constr = 0;
    int cnt     = 0;
    for (;;){
        skipWhitespace(in);
        if (*in == EOF) break;
        else if (*in == 'p'){
            if (eagerMatch(in, "p cnf")){ // || eagerMatch(in, "p cnf+")){
                if(*in == '+') ++in;
                vars    = parseInt(in);
                while (vars > S.nVars()) S.newVar();  // Caution: this will break on some
                                                      // malformed instances (with incorrect
                                                      // headers)
                constr = parseInt(in);
                // SATRACE'06 hack
                // if (clauses > 4000000)
                //     S.eliminate(true);
            }
            else {
                printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
            }
        } else if (*in == 'c' || *in == 'p')
            skipLine(in);
        else {
            cnt++;
            if(readConstr(in, S, lits, bound)) S.addAtMost_(lits,bound);
            else S.addClause_(lits);
        }
    }
    if (vars != S.nVars())
        fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of variables.\n");
    if (cnt  != constr)
        fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of constraints.\n");
}

// Inserts problem into solver.
//
template<class Solver>
static void parse_DIMACS(gzFile input_stream, Solver& S) {
    StreamBuffer in(input_stream);
    parse_DIMACS_main(in, S); }

//=================================================================================================
}

#endif
