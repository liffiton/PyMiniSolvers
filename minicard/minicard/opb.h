/****************************************************************************************[opb.h]
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

#ifndef Minisat_opb_h
#define Minisat_opb_h

#include <stdio.h>

#include "utils/ParseUtils.h"
#include "minicard/SolverTypes.h"

namespace Minisat {

//=================================================================================================
// OPB Parser:

template<class B, class Solver>
static void readConstr(B& in, Solver& S) {
    vec<Lit> lits;
    int weight = 0;
    int boundoffset = 0;

    for (;;){
        skipWhitespace(in);
        if (eagerMatch(in, "<=")) {
            int bound = parseInt(in);
            bound += boundoffset;
            S.addAtMost_(lits, bound);
            return;
        }
        else if (eagerMatch(in, ">=")) {
            int bound = parseInt(in);
            bound += boundoffset;
            if (bound == 1) {
                S.addClause_(lits);
            }
            else {
                bound = lits.size() - bound;
                for (int i=0;i<lits.size();i++) {
                    lits[i] = ~lits[i];
                }
                S.addAtMost_(lits, bound);
            }
            return;
        }
        else if (eagerMatch(in, "=")) {
            int bound = parseInt(in);
            bound += boundoffset;
            S.addAtMost_(lits, bound);
            bound = lits.size() - bound;
            for (int i=0;i<lits.size();i++) {
                lits[i] = ~lits[i];
            }
            S.addAtMost_(lits, bound);
            return;
        }
        else if (eagerMatch(in, "x")) {
            int parsed_lit = parseInt(in);
            int var = abs(parsed_lit)-1; 
            while (var >= S.nVars()) S.newVar();
            assert(weight != 0);
            lits.push((weight < 0) ? ~mkLit(var) : mkLit(var));
        }
        else {
            weight = parseInt(in);
            if (weight < 0) {
                // -k * x_i ... >= b  --->  ~x_i ... >= b+k
                boundoffset -= weight;
            }
        }
    }
}

template<class B, class Solver>
static void parse_OPB_main(B& in, Solver& S) {
    for (;;){
        skipWhitespace(in);
        if (*in == EOF) {
            break;
        }
        else if (*in == '*' || *in == ';') {
            // comment or end of a line
            skipLine(in);
        }
        else {
            readConstr(in, S);
            //S.addAtMost_(lits,bound);
        }
    }
}

// Inserts problem into solver.
//
template<class Solver>
static void parse_OPB(gzFile input_stream, Solver& S) {
    StreamBuffer in(input_stream);
    parse_OPB_main(in, S); }

//=================================================================================================
}

#endif
