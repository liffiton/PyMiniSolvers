"""A Python API for the MiniSat_ and MiniCard_ constraint solvers.

.. _MiniSat: http://minisat.se/
.. _MiniCard: http://git.io/minicard

Classes:
  Solver
    An abstract base class for the other classes.
  SubsetMixin
    A mixin class adding 'subset' functionality to Solver subclasses.

  :class:`MinisatSolver`
    Solve CNF instances using MiniSat.
  :class:`MinicardSolver`
    Solve CNF+ (CNF plus cardinality constraints) using MiniCard.

  :class:`MinisatSubsetSolver`
    Solve arbitrary subsets of CNF instances and find SAT subsets / UNSAT cores.
  :class:`MinicardSubsetSolver`
    Solve arbitrary subsets of CNF+ instances and find SAT subsets / UNSAT cores.
"""

import array
import os
import ctypes
from abc import ABCMeta, abstractmethod
from ctypes import c_void_p, c_ubyte, c_bool, c_int


class Solver(object):
    """The Solver class is an abstract base class for MiniSat and
    MiniCard solver classes.  It provides the basic methods that both
    contain, closely following the methods in MiniSat and MiniCard's
    Solver class.

    Solver should not be instantiated directly.  Instead, use its
    subclasses MinisatSolver, MinicardSolver, MinisatSubsetSolver, or
    MinicardSubsetSolver (see below).
    """

    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, libfilename):
        self._setup_lib(libfilename)
        self.s = self.lib.Solver_new()

    def _setup_lib(self, libfilename):
        """Load the minisat library with ctypes and create a Solver
           object.  Correct return types (if not int as assumed by
           ctypes) and set argtypes for functions from the minisat
           library.
        """
        dirname = os.path.dirname(os.path.abspath(__file__))
        libfile = dirname + '/' + libfilename
        if not os.path.exists(libfile):
            raise IOError("Specified library file not found.  Did you run 'make' to build the solver libraries?\nFile not found: %s" % libfile)

        self.lib = ctypes.cdll.LoadLibrary(dirname+'/'+libfilename)

        l = self.lib

        l.Solver_new.restype = c_void_p
        l.Solver_new.argtypes = []
        l.Solver_delete.argtypes = [c_void_p]

        l.nVars.argtypes = [c_void_p]
        l.nClauses.argtypes = [c_void_p]
        l.setPhaseSaving.argtypes = [c_void_p, c_int]
        l.setRndPol.argtypes = [c_void_p, c_bool]

        l.newVar.argtypes = [c_void_p, c_ubyte, c_bool]

        l.addClause.restype = c_bool
        l.addClause.argtypes = [c_void_p, c_int, c_void_p]
        l.addUnit.restype = c_bool
        l.addUnit.argtypes = [c_void_p, c_int]

        l.solve.restype = c_bool
        l.solve.argtypes = [c_void_p]
        l.solve_assumptions.restype = c_bool
        l.solve_assumptions.argtypes = [c_void_p, c_int, c_void_p]
        l.simplify.restype = c_bool
        l.simplify.argtypes = [c_void_p]

        l.unsatCore.argtypes = [c_void_p, c_int, c_void_p]
        l.modelValue.argtypes = [c_void_p, c_int]
        l.fillModel.argtypes = [c_void_p, c_void_p, c_int, c_int]
        l.getModelTrues.restype = c_int
        l.getModelTrues.argtypes = [c_void_p, c_void_p, c_int, c_int]

        l.getImplies.argtypes = [c_void_p, c_void_p]
        l.getImplies.restype = c_int

    def __del__(self):
        """Delete the Solver object"""
        self.lib.Solver_delete(self.s)

    @staticmethod
    def _to_intptr(a):
        """Helper function to get a ctypes POINTER(c_int) for an array"""
        addr, size = a.buffer_info()
        return ctypes.cast(addr, ctypes.POINTER(c_int)), size

    def new_var(self, polarity=None, dvar=True):
        """Create a new variable in the solver.

        Args:
            polarity (bool):
              The default polarity for this variable.  True = variable's
              default is True, etc.  Note that this is the reverse of the 'user
              polarity' in MiniSat, where True indicates the *sign* is True,
              hence the default value is False.
            dvar (bool):
              Whether this variable will be used as a decision variable.

        Returns:
            The new variable's index (0-based counting).
        """

        if polarity is None:
            pol_int = 2
        elif polarity is True:
            pol_int = 1
        elif polarity is False:
            pol_int = 0
        return self.lib.newVar(self.s, pol_int, dvar)

    def nvars(self):
        return self.lib.nVars(self.s)

    def nclauses(self):
        return self.lib.nClauses(self.s)

    def set_phase_saving(self, ps):
        '''Set the level of phase saving (0=none, 1=limited, 2=full (default)).'''
        self.lib.setPhaseSaving(self.s, ps)

    def set_rnd_pol(self, val):
        '''Set whether random polarities are used for decisions (overridden if vars are created with a user polarity other than None)'''
        self.lib.setRndPol(self.s, val)

    def add_clause(self, lits):
        """Add a clause to the solver.

        Args:
            lits:
              A list of literals as integers.  Each integer specifies a
              variable with *1*-based counting and a sign via the sign of the
              integer.  Ex.: [-1, 2, -3] is (!x0 + x1 + !x2)

        Returns:
            A boolean value returned from MiniSat's ``addClause()`` function,
            indicating success (True) or conflict (False).
        """
        if not all(abs(x) <= self.nvars() for x in lits):
            raise Exception("Not all variables in %s are created yet.  Call new_var() first." % lits)
        if len(lits) > 1:
            a = array.array('i', lits)
            a_ptr, size = self._to_intptr(a)
            return self.lib.addClause(self.s, size, a_ptr)
        elif len(lits) == 1:
            return self.lib.addUnit(self.s, lits[0])
        else:
            return self.lib.addClause(self.s, 0, None)

    def solve(self, assumptions=None):
        """Solve the current set of clauses, optionally with a set of assumptions.

        Args:
            assumptions:
              An iterable returning literals as integers, specified as in
              ``add_clause()``.

        Returns:
            True if the clauses (and assumptions) are satisfiable, False otherwise.
        """
        if assumptions is None:
            return self.lib.solve(self.s)
        else:
            a = array.array('i', assumptions)
            a_ptr, size = self._to_intptr(a)
            return self.lib.solve_assumptions(self.s, size, a_ptr)

    def simplify(self):
        return self.lib.simplify(self.s)

    def get_model(self, start=0, end=-1):
        """Get the current model from the solver, optionally retrieving only a slice.

        Args:
            start, end (int):
              Optional start and end indices, interpreted as in ``range()``.

        Returns:
            An array of booleans indexed to each variable (from 0).  If a start
            index was given, the returned list starts at that index (i.e.,
            ``get_model(10)[0]`` is index 10 from the solver's model.
        """
        if end == -1:
            end = self.nvars()
        a = array.array('i', [-1] * (end-start))
        a_ptr, size = self._to_intptr(a)
        self.lib.fillModel(self.s, a_ptr, start, end)
        return a

    def get_model_trues(self, start=0, end=-1):
        """Get variables assigned true in the current model from the solver.

        Args:
            start, end (int):
              Optional start and end indices, interpreted as in ``range()``.

        Returns:
            An array of true variables in the solver's current model.  If a
            start index was given, the variables are indexed from that value.
        """
        if end == -1:
            end = self.nvars()
        a = array.array('i', [-1] * (end-start))
        a_ptr, size = self._to_intptr(a)
        count = self.lib.getModelTrues(self.s, a_ptr, start, end)
        # reduce the array down to just the valid indexes
        return a[:count]

    def model_value(self, i):
        return self.lib.modelValue(self.s, i)

    def implies(self):
        """Get literals known to be implied by the current formula.
           (I.e., all assignments made at level 0)

        Returns:
           An array of literals.
        """
        a = array.array('i', [-1] * self.nvars())
        a_ptr, size = self._to_intptr(a)
        count = self.lib.getImplies(self.s, a_ptr)
        # reduce the array down to just the valid indexes
        return a[:count]


class SubsetMixin(object):
    """A mixin for any Solver class that lets it reason about subsets of a clause set."""
    _origvars = None
    _relvars = None

    def set_varcounts(self, vars, constraints):
        """Record how many of the solver's variables and clauses are
           "original," as opposed to clause-selector variables, etc.
        """
        self._origvars = vars
        self._relvars = constraints

    def add_clause_instrumented(self, lits, index):
        """Add a clause with a relaxation variable (the rel.var. is
           based on the index, which is assumed to be 0-based).
        """
        if self._origvars is None:
            raise Exception("SubsetSolver.set_varcounts() must be called before .add_clause_instrumented()")
        instrumented_clause = [-(self._origvars+1+index)] + lits
        self.add_clause(instrumented_clause)

    def solve_subset(self, subset):
        if self._origvars is None:
            raise Exception("SubsetSolver.set_varcounts() must be called before .solve_subset()")
        # convert clause indices to clause-selector variable indices
        a = array.array('i', (i+self._origvars+1 for i in subset))
        a_ptr, size = self._to_intptr(a)
        return self.lib.solve_assumptions(self.s, size, a_ptr)

    def unsat_core(self):
        a = array.array('i', [-1] * self.nclauses())
        a_ptr, size = self._to_intptr(a)
        length = self.lib.unsatCore(self.s, self._origvars, a_ptr)
        # reduce the array down to just the valid indexes
        return a[:length]

    def sat_subset(self):
        return self.get_model_trues(start=self._origvars, end=self._origvars+self._relvars)


class MinisatSolver(Solver):
    """A Python analog to MiniSat's Solver class.

    >>> S = MinisatSolver()

    Create variables using ``new_var()``.  Add clauses as list of literals with
    ``add_clause()``, analogous to MiniSat's ``add_clause()``.  Literals are
    specified as integers, with the magnitude indicating the variable index
    (with 1-based counting) and the sign indicating True/False.  For example,
    to add clauses (x0), (!x1), (!x0 + x1 + !x2), and (x2 + x3):

    >>> for i in range(4):
    ...     S.new_var()  # doctest: +ELLIPSIS
    0
    1
    ...
    >>> for clause in [1], [-2], [-1, 2, -3], [3, 4]:
    ...     S.add_clause(clause)  # doctest: +ELLIPSIS
    True
    True
    ...

    The ``solve()`` method returns True or False just like MiniSat's.

    >>> S.solve()
    True

    Models are returned as arrays of Booleans, indexed by var.
    So the following represents x0=True, x1=False, x2=False, x3=True.

    >>> list(S.get_model())
    [1, 0, 0, 1]

    The ``add_clause()`` method may return False if a conflict is detected
    when adding the clause, even without search.

    >>> S.add_clause([-4])
    False
    >>> S.solve()
    False
    """
    def __init__(self):
        super(MinisatSolver, self).__init__("libminisat.so")


class MinicardSolver(Solver):
    """A Python analog to MiniCard's Solver class.

    >>> S = MinicardSolver()
    >>> for i in range(4):  tmp=S.new_var()

    Add clauses (x0), (!x1), and (x2 + x3)
    >>> for clause in [1], [-2], [3, 4]:
    ...    S.add_clause(clause)
    True
    True
    True

    AtMost({x0, !x1, x2}, 2)
    >>> S.add_atmost([1,-2,3], 2)
    True

    >>> S.solve()
    True

    Models are returned as arrays of Booleans, indexed by var.
    So the following represents x0=True, x1=False, x2=False, x3=True.
    >>> list(S.get_model())
    [1, 0, 0, 1]
    """
    def __init__(self):
        super(MinicardSolver, self).__init__("libminicard.so")

    def _setup_lib(self, libfilename):
        """Correct return types (if not int as assumed by ctypes) and set argtypes for
           functions from the minicard library.
        """
        super(MinicardSolver, self)._setup_lib(libfilename)

        # additional function for minicard
        l = self.lib
        l.addAtMost.restype = c_bool
        l.addAtMost.argtypes = [c_void_p, c_int, c_void_p, c_int]

    def add_atmost(self, lits, k):
        """Add an AtMost constraint to the solver.

        Args:
            lits:
              A list of literals as integers.  Each integer specifies a
              variable with **1**-based counting and a sign via the sign of
              the integer.  Ex.: [-1, 2, -3] is {!x0, x1, !x2}
            k (int):
              The [upper] bound to place on these literals.

        Returns:
            A boolean value returned from MiniCard's ``addAtMost()``
            function, indicating success (True) or conflict (False).
        """
        if not all(abs(x) <= self.nvars() for x in lits):
            raise Exception("Not all variables in %s are created yet.  Call new_var() first." % lits)

        if len(lits) > 1:
            a = array.array('i', lits)
            a_ptr, size = self._to_intptr(a)
            return self.lib.addAtMost(self.s, size, a_ptr, k)
        else:
            return self.lib.addAtMost(self.s, 0, None, 0)


class MinisatSubsetSolver(SubsetMixin, MinisatSolver):
    """A class for reasoning about subsets of constraints within MiniSat.

    >>> S = MinisatSubsetSolver()

    It must be told explicitlyhow many of its variables are "real" and how
    many are relaxation variables for constraints.
    >>> S.set_varcounts(vars = 4, constraints = 5)

    >>> for i in range(4+5):  tmp=S.new_var()
    >>> for i, clause in enumerate([[1], [-2], [-1, 2, 3], [-3], [-1]]):
    ...    S.add_clause_instrumented(clause, i)

    Any subset of the constraints can be tested for satisfiability.
    >>> S.solve_subset([0,1,2])
    True
    >>> S.solve_subset([0,1,2,3])
    False

    If a subset is found to be unsatisfiable, an UNSAT core can be found.
    Cores are returned as array objects.
    >>> core = S.unsat_core()
    >>> sorted(core)
    [0, 1, 2, 3]
    """

    pass


class MinicardSubsetSolver(SubsetMixin, MinicardSolver):
    """A class for reasoning about subsets of constraints within
    MiniCard.

    >>> S = MinicardSubsetSolver()
    >>> S.set_varcounts(vars = 4, constraints = 4)
    >>> for i in range(4+4):  tmp=S.new_var()
    >>> for i, clause in enumerate([[1], [-2], [3], [4]]):
    ...    S.add_clause_instrumented(clause, i)

    AtMost constraints cannot be instrumented -- they are simply hard
    constraints.
    >>> S.add_atmost([1,-2,3], 2)
    True

    Any subset of the constraints can be tested for satisfiability.
    >>> S.solve_subset([0,1])
    True
    >>> S.solve_subset([0,1,2,3])
    False

    If a subset is found to be unsatisfiable, an UNSAT core can be
    found.  Cores are returned as array objects.  Hard constraints are
    not returned in the core, but the core is found with respect to those
    constraints as well.
    >>> core = S.unsat_core()
    >>> sorted(core)
    [0, 1, 2]
    """

    pass
