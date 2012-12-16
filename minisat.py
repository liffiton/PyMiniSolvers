from ctypes import *
import os

class Solver:
    def __init__(self):
        dir = os.path.dirname(os.path.abspath(__file__))
        self.lib = cdll.LoadLibrary(dir+'/libminisat.so')
        self.setup_lib()
        self.s = self.lib.Solver_new()
        self.model = []
        self.model_stale = True

    def setup_lib(self):
        l = self.lib
        # correct return types (if not the assumed int) and set argtypes
        l.Solver_new.restype = c_void_p
        l.Solver_new.argtypes = []
        l.Solver_delete.argtypes = [c_void_p]

        l.nVars.argtypes = [c_void_p]
        l.nClauses.argtypes = [c_void_p]

        l.newVar.argtypes = [c_void_p, c_ubyte]

        l.addClause.restype = c_bool
        l.addClause.argtypes = [c_void_p, c_int, c_void_p]
        l.addUnit.restype = c_bool
        l.addUnit.argtypes = [c_void_p, c_int]

        l.solve.restype = c_bool
        l.solve.argtypes = [c_void_p]
        l.solve_assumptions.restype = c_bool
        l.solve_assumptions.argtypes = [c_void_p, c_int, c_void_p]
        l.solve_subset.restype = c_bool
        l.solve_subset.argtypes = [c_void_p, c_int, c_void_p]
        l.simplify.restype = c_bool
        l.simplify.argtypes = [c_void_p]

        l.unsatCore.argtypes = [c_void_p, c_void_p]
        l.modelValue.argtypes = [c_void_p, c_int]
        l.fillModel.argtypes = [c_void_p, c_void_p]

    def __del__(self):
        self.lib.Solver_delete(self.s)

    def new_var(self, polarity=None):
        if polarity is None:
            pol_int = 2
        elif polarity == True:
            pol_int = 1
        elif polarity == False:
            pol_int = 0
        self.lib.newVar(self.s, pol_int)

    def nvars(self):
        return self.lib.nVars(self.s)

    def nclauses(self):
        return self.lib.nClauses(self.s)

    def add_clause(self, lits):
        if len(lits) > 1:
            array = (c_int * len(lits))(*lits)
            self.lib.addClause(self.s, len(lits), array)
        elif len(lits) == 1:
            self.lib.addUnit(self.s, lits[0])
        else:
            self.lib.addClause(self.s, 0, None)

    def solve(self, assumptions=None):
        self.model_stale = True
        if assumptions is None:
            return self.lib.solve(self.s)
        else:
            array = (c_int * len(assumptions))(*assumptions)
            return self.lib.solve_assumptions(self.s, len(assumptions), array)

    def solve_subset(self, subset):
        self.model_stale = True
        array = (c_int * len(subset))(*subset)
        return self.lib.solve_subset(self.s, len(subset), array)

    def simplify(self):
        return self.lib.simplify(self.s)

    def update_model(self):
        array = (c_int * self.nvars())()
        self.lib.fillModel(self.s, array)
        self.model = array[:]
        self.model_stale = False

    def get_model(self):
        if self.model_stale:
            self.update_model()
        return self.model

    def model_value(self, i):
        if self.model_stale:
            self.update_model()
        return self.model[i-1]
        # slower:
        #return self.lib.modelValue(self.s, i)

    def unsat_core(self):
        array = (c_int * self.nclauses())()
        length = self.lib.unsatCore(self.s, array)
        return array[:length]

