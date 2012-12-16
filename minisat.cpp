#include "minisat/core/Solver.h"

using namespace Minisat;

Lit itoLit(int i) {
    bool sign = i < 0;
    int var = (sign) ? -i-1 : i-1; // 0-based variable numbering
    return (sign) ? ~mkLit(var) : mkLit(var);
}

extern "C" {
    Solver* Solver_new() { return new Solver(); }
    void Solver_delete(Solver* s) { delete s; }

    int nVars(Solver* s) { return s->nVars(); }
    int nClauses(Solver* s) { return s->nClauses(); }

    void newVar(Solver* s, uint8_t polarity) { s->newVar(lbool(polarity)); }  // 0=False, 1=True, 2=Undef
    bool addClause(Solver* s, int len, int* lits) {
        vec<Lit> clause;
        for (int i = 0 ; i < len ; i++) {
            clause.push( itoLit(lits[i]) );
        }
        return s->addClause(clause);
    }
    bool addUnit(Solver* s, int lit) {
        return s->addClause(itoLit(lit));
    }

    bool solve(Solver* s) { return s->solve(); }
    bool solve_assumptions(Solver* s, int len, int* lits) {
        vec<Lit> assumptions;
        for (int i = 0 ; i < len ; i++) {
            assumptions.push( itoLit(lits[i]) );
        }
        return s->solve(assumptions);
    }
    bool solve_subset(Solver* s, int len, int* subset) {
        vec<Lit> assumptions;
        int nv = s->nVars() - s->nClauses();
        for (int i = 0 ; i < len ; i++) {
            assumptions.push( itoLit(subset[i] + nv + 1) );
        }
        return s->solve(assumptions);
    }

    bool simplify(Solver* s) {
        return s->simplify();
    }

    // This is fairly slow to call from Python.
    // It is better to copy the whole model over with fillModel()
    // if you will be looking at most or all values, anyway.
    int modelValue(Solver* s, int i) { return s->modelValue(i-1) != l_False; }

    void fillModel(Solver* s, int* model) {
        for (int i = 0 ; i < s->nVars() ; i++) {
            model[i] = s->modelValue(i) != l_False;
        }
    }

    // returns a core with 0-based counting
    // (i.e., first clause is 0, etc.)
    int unsatCore(Solver* s, int* core) {
        int nv = s->nVars() - s->nClauses();
        for (int i = 0 ; i < s->conflict.size() ; i++) {
            core[i] = var(s->conflict[i]) - nv;
        }
        return s->conflict.size();
    }
}
