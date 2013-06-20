import array
import os
import ctypes
from ctypes import c_void_p, c_ubyte, c_bool, c_int

class Solver(object):
    """The Solver class is meant as a fairly direct analog of Minisat/Minicard's Solver class."""

    def __init__(self, libfilename):
        """Load the minisat library with ctypes and create a Solver object."""
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.lib = ctypes.cdll.LoadLibrary(dirname+'/'+libfilename)
        self._setup_lib()
        self.s = self.lib.Solver_new()

    def _setup_lib(self):
        """Correct return types (if not int as assumed by ctypes) and set argtypes for
           functions from the minisat library.
        """
        l = self.lib

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
        l.solve_subset.argtypes = [c_void_p, c_int, c_int, c_void_p]
        l.simplify.restype = c_bool
        l.simplify.argtypes = [c_void_p]

        l.unsatCore.argtypes = [c_void_p, c_int, c_void_p]
        l.modelValue.argtypes = [c_void_p, c_int]
        l.fillModel.argtypes = [c_void_p, c_void_p, c_int, c_int]
        l.getModelTrues.restype = c_int
        l.getModelTrues.argtypes = [c_void_p, c_void_p, c_int, c_int]

    def __del__(self):
        """Delete the Solver object"""
        self.lib.Solver_delete(self.s)

    @staticmethod
    def _to_intptr(a):
        """Helper function to get a ctypes POINTER(c_int) for an array"""
        addr, size = a.buffer_info()
        return ctypes.cast(addr, ctypes.POINTER(c_int)) , size

    def new_var(self, polarity=None):
        """Create a new variable in the solver.
        
        Args:
            polarity: A boolean specifying the default polarity for this
                variable.  True = variable's default is True, etc.  Note
                that this is the reverse of the 'user polarity' in minisat,
                where True indicates the *sign* is True, hence the default
                value is False.

        Returns:
            The new variable's number (0-based counting).
        """
              
        if polarity is None:
            pol_int = 2
        elif polarity == True:
            pol_int = 1
        elif polarity == False:
            pol_int = 0
        return self.lib.newVar(self.s, pol_int)

    def nvars(self):
        return self.lib.nVars(self.s)

    def nclauses(self):
        return self.lib.nClauses(self.s)

    def add_clause(self, lits):
        """Add a clause to the solver.

        Args:
            lits: A list of literals as integers.  Each integer specifies a
                variable with *1*-based counting and a sign via the sign of
                the integer.  Ex.: [-1, 2, -3] is [!x0, x1, !x2]

        Returns:
            A boolean value returned from Minisat's addClause() function,
            indicating success (True) or conflict (False).
        """
        if not all(abs(x) <= self.nvars() for x in lits):
            raise Exception("Not all variables in %s are created yet.  Call new_var() first." % lits)
        if len(lits) > 1:
            a = array.array('i',lits)
            a_ptr, size = self._to_intptr(a)
            self.lib.addClause(self.s, size, a_ptr)
        elif len(lits) == 1:
            self.lib.addUnit(self.s, lits[0])
        else:
            self.lib.addClause(self.s, 0, None)

    def solve(self, assumptions=None):
        """Solve the current set of clauses, optionally with a set of assumptions.

        Args:
            assumptions: An iterable returning literals as integers, specified as in add_clause().

        Returns:
            True if the clauses (and assumptions) are satisfiable, False otherwise.
        """
        if assumptions is None:
            return self.lib.solve(self.s)
        else:
            a = array.array('i',assumptions)
            a_ptr, size = self._to_intptr(a)
            return self.lib.solve_assumptions(self.s, size, a_ptr)

    def simplify(self):
        return self.lib.simplify(self.s)

    def get_model(self, start=0, end=-1):
        """Get the current model from the solver, optionally retrieving only a slice.

        Args:
            start, end:  The optional start and end indices, interpreted as in range().

        Returns:
            A list of booleans indexed to each variable (from 0).  If a start index was
            given, the returned list starts at that index (i.e., get_model(10)[0] is index
            10 from the solver's model.
        """
        if end == -1:
            end = self.nvars()
        a = (c_int * (end-start))()
        self.lib.fillModel(self.s, a, start, end)
        return a[:]

    def get_model_trues(self, start=0, end=-1):
        """Get variables assigned true in the current model from the solver.

        Args:
            start, end:  The optional start and end indices, interpreted as in range().

        Returns:
            A list of true variables in the solver's current model.  If a start
            index was given, the variables are indexed from that value.
        """
        if end == -1:
            end = self.nvars()
        a = (c_int * (end-start))()
        count = self.lib.getModelTrues(self.s, a, start, end)
        return a[:count]

    def model_value(self, i):
        return self.lib.modelValue(self.s, i)

class SubsetMixin(object):
    """A mixin for any Solver class that lets it reason about subsets of a clause set."""
    origvars = None
    origclauses = None
    n = 0

    def set_orig(self, o_vars, o_clauses):
        """Record how many of the solver's variables and clauses are "original,"
            as opposed to clause-selector variables, etc.
        """
        self.origvars = o_vars
        self.origclauses = o_clauses

    def add_clause(self, lits):
        if self.origvars is None:
            raise Exception("SubsetSolver.set_orig() must be called before .add_clause()")
        instrumented_clause = [-(self.origvars+self.n+1)] + lits
        super(SubsetMixin, self).add_clause(instrumented_clause)
        self.n += 1

    def solve_subset(self, subset):
        a = array.array('i', subset)
        a_ptr, size = self._to_intptr(a)
        return self.lib.solve_subset(self.s, self.origvars, size, a_ptr)

    def unsat_core(self):
        a = (c_int * self.nclauses())()
        length = self.lib.unsatCore(self.s, self.origvars, a)
        return a[:length]

    def sat_subset(self):
        return self.get_model_trues(start = self.origvars, end = self.origvars+self.origclauses)

class MinisatSolver(Solver):
    def __init__(self):
        super(MinisatSolver, self).__init__("libminisat.so")

class MinicardSolver(Solver):
    def __init__(self):
        super(MinicardSolver, self).__init__("libminicard.so")

    def _setup_lib(self):
        """Correct return types (if not int as assumed by ctypes) and set argtypes for
           functions from the minicard library.
        """
        super(MinicardSolver, self)._setup_lib()

        # additional function for minicard
        l = self.lib
        l.addAtMost.restype = c_bool
        l.addAtMost.argtypes = [c_void_p, c_int, c_void_p, c_int]

    def add_atmost(self, lits, k):
        """Add an AtMost constraint to the solver.

        Args:
            lits: A list of literals as integers.  Each integer specifies a
                variable with *1*-based counting and a sign via the sign of
                the integer.  Ex.: [-1, 2, -3] is [!x0, x1, !x2]
            k: The [upper] bound to place on these literals

        Returns:
            A boolean value returned from MiniCard's addAtMost() function,
            indicating success (True) or conflict (False).
        """
        if not all(abs(x) <= self.nvars() for x in lits):
            raise Exception("Not all variables in %s are created yet.  Call new_var() first." % lits)
        if len(lits) > 1:
            array = (c_int * len(lits))(*lits)
            self.lib.addAtMost(self.s, len(lits), array, k)
        else:
            self.lib.addAtMost(self.s, 0, None, 0)

class MinisatSubsetSolver(SubsetMixin, MinisatSolver):
    pass

class MinicardSubsetSolver(SubsetMixin, MinicardSolver):
    pass

