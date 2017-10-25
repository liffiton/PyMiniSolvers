"""A Python API for the MiniSat_ and MiniCard_ constraint solvers.

.. _MiniSat: http://minisat.se/
.. _MiniCard: http://git.io/minicard

Classes:
  `MinisatSolver`
    Solve CNF instances using MiniSat.
  `MinicardSolver`
    Solve CNF+ (CNF plus cardinality constraints) using MiniCard.

  `MinisatSubsetSolver`
    Solve arbitrary subsets of CNF instances and find SAT subsets / UNSAT cores.
  `MinicardSubsetSolver`
    Solve arbitrary subsets of CNF+ instances and find SAT subsets / UNSAT cores.

  Solver
    An abstract base class for the other classes.
  SubsetMixin
    A mixin class adding 'subset' functionality to Solver subclasses.
"""

import array
import os
import ctypes  # type: ignore
from abc import ABCMeta, abstractmethod
from ctypes import c_void_p, c_ubyte, c_bool, c_int, c_int64, c_double  # type: ignore

try:
    import typing  # noqa: for mypy-lang type-checking
    from typing import Iterable, Sequence, Tuple  # noqa: for mypy-lang type-checking
except ImportError:
    # not needed at runtime, so no error
    pass


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
    def __init__(self, libfilename):  # type: (str) -> None
        self._setup_lib(libfilename)
        self.s = self.lib.Solver_new()

    def _setup_lib(self, libfilename):  # type: (str) -> None
        """Load the minisat library with ctypes and create a Solver
           object.  Correct return types (if not int as assumed by
           ctypes) and set argtypes for functions from the minisat
           library.
        """
        dirname = os.path.dirname(os.path.abspath(__file__))
        libfile = os.path.join(dirname, libfilename)
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
        l.setRndInitAct.argtypes = [c_void_p, c_bool]
        l.setRndSeed.argtypes = [c_void_p, c_double]

        l.newVar.argtypes = [c_void_p, c_ubyte, c_bool]

        l.addClause.restype = c_bool
        l.addClause.argtypes = [c_void_p, c_int, c_void_p]
        l.addUnit.restype = c_bool
        l.addUnit.argtypes = [c_void_p, c_int]

        l.solve.restype = c_bool
        l.solve.argtypes = [c_void_p]
        l.solve_assumptions.restype = c_bool
        l.solve_assumptions.argtypes = [c_void_p, c_int, c_void_p]
        l.check_complete.restype = c_bool
        l.check_complete.argtypes = [c_void_p, c_int, c_void_p, c_bool]
        l.simplify.restype = c_bool
        l.simplify.argtypes = [c_void_p]

        l.conflictSize.argtypes = [c_void_p]
        l.conflictSize.restype = c_int
        l.unsatCore.argtypes = [c_void_p, c_int, c_void_p, c_int]
        l.unsatCore.restype = c_int
        l.modelValue.argtypes = [c_void_p, c_int]
        l.modelValue.restype = c_int
        l.fillModel.argtypes = [c_void_p, c_void_p, c_int, c_int]
        l.getModelTrues.argtypes = [c_void_p, c_void_p, c_int, c_int, c_int]
        l.getModelTrues.restype = c_int

        l.getImplies.argtypes = [c_void_p, c_void_p]
        l.getImplies.restype = c_int
        l.getImplies_assumptions.argtypes = [c_void_p, c_void_p, c_void_p, c_int]
        l.getImplies_assumptions.restype = c_int

        l.get_solves.argtypes = [c_void_p]
        l.get_solves.restype = c_int64
        l.get_starts.argtypes = [c_void_p]
        l.get_starts.restype = c_int64
        l.get_decisions.argtypes = [c_void_p]
        l.get_decisions.restype = c_int64
        l.get_rnd_decisions.argtypes = [c_void_p]
        l.get_rnd_decisions.restype = c_int64
        l.get_propagations.argtypes = [c_void_p]
        l.get_propagations.restype = c_int64
        l.get_conflicts.argtypes = [c_void_p]
        l.get_conflicts.restype = c_int64

    def __del__(self):  # type: () -> None
        """Delete the Solver object"""
        self.lib.Solver_delete(self.s)

    @staticmethod
    def _to_intptr(a):  # type: (array.array) -> Tuple[int, int]
        """Helper function to get a ctypes POINTER(c_int) for an array"""
        addr, size = a.buffer_info()
        return ctypes.cast(addr, ctypes.POINTER(c_int)), size

    @staticmethod
    def _get_array(seq):  # type: (Iterable[int]) -> array.array
        """Helper function to turn any iterable into an array (unless it already is one)"""
        if isinstance(seq, array.array):
            return seq
        else:
            return array.array('i', seq)

    def new_var(self, polarity=None, dvar=True):  # type: (bool, bool) -> int
        """Create a new variable in the solver.

        Args:
            polarity (bool):
              The default polarity for this variable.  True = variable's
              default is True, etc.  Note that this is the reverse of the 'user
              polarity' in MiniSat, where True indicates the *sign* is True.
              The default, None, creates the variable using Minisat's default,
              which assigns a variable False at first, but then may change that
              based on the phase-saving setting.
            dvar (bool):
              Whether this variable will be used as a decision variable.

        Returns:
            The new variable's index (0-based counting).
        """

        if polarity is None:
            pol_int = 2  # lbool l_Undef
        elif polarity is True:
            pol_int = 1  # lbool l_False (hence, the *sign* is false, so the literal is true)
        elif polarity is False:
            pol_int = 0  # lbool l_True (hence the literal is false)
        return self.lib.newVar(self.s, pol_int, dvar)

    def nvars(self):  # type: () -> int
        '''Get the number of variables created in the solver.'''
        return self.lib.nVars(self.s)

    def nclauses(self):  # type: () -> int
        '''Get the number of clauses or constraints added to the solver.'''
        return self.lib.nClauses(self.s)

    def set_phase_saving(self, ps):  # type: (int) -> None
        '''Set the level of phase saving (0=none, 1=limited, 2=full (default)).'''
        self.lib.setPhaseSaving(self.s, ps)

    def set_rnd_pol(self, val):  # type: (bool) -> None
        '''Set whether random polarities are used for decisions (overridden if vars are created with a user polarity other than None)'''
        self.lib.setRndPol(self.s, val)

    def set_rnd_init_act(self, val):  # type: (bool) -> None
        '''Set whether variables are intialized with a random initial activity.
           (default: False)'''
        self.lib.setRndInitAct(self.s, val)

    def set_rnd_seed(self, seed):  # type: (float) -> None
        '''Set the solver's random seed to the given double value.  Cannot be 0.0.'''
        assert(seed != 0.0)
        self.lib.setRndSeed(self.s, seed)

    def add_clause(self, lits):  # type: (Sequence[int]) -> bool
        """Add a clause to the solver.

        Args:
            lits:
              A sequence of literals as integers.  Each integer specifies a
              variable with *1*-based counting and a sign via the sign of the
              integer.  Ex.: [-1, 2, -3] is (!x0 + x1 + !x2)

        Returns:
            A boolean value returned from MiniSat's ``addClause()`` function,
            indicating success (True) or conflict (False).
        """
        if not all(abs(x) <= self.nvars() for x in lits):
            raise Exception("Not all variables in %s are created yet.  Call new_var() first." % lits)
        if len(lits) > 1:
            a = self._get_array(lits)
            a_ptr, size = self._to_intptr(a)
            return self.lib.addClause(self.s, size, a_ptr)
        elif len(lits) == 1:
            (lit,) = lits   # extract one item whether list or set
            return self.lib.addUnit(self.s, lit)
        else:
            return self.lib.addClause(self.s, 0, None)

    def check_complete(self, positive_lits=None, negative_lits=None):  # type: (Sequence[int], Sequence[int]) -> bool
        """Check whether a given complete assignment satisfies the current set
        of clauses.  For efficiency, it may be given just the positive literals
        or just the negative literals.

        Args:
            positive_lits, negative_lits:
              Optional sequences (exactly one must be specified) containing
              literals as integers, specified as in `add_clause()`.  If
              positive literals are given, the assignment will be completed
              assuming all other variables are negative, and vice-versa if
              negative literals are given.

        Returns:
            True if the assignment satisfies the current clauses, False otherwise.
        """
        if positive_lits is not None:
            a = self._get_array(positive_lits)
            a_ptr, size = self._to_intptr(a)
            return self.lib.check_complete(self.s, size, a_ptr, True)
        elif negative_lits is not None:
            a = self._get_array(negative_lits)
            a_ptr, size = self._to_intptr(a)
            return self.lib.check_complete(self.s, size, a_ptr, False)
        else:
            raise Exception("Either positive_lits or negative_lits must be specified in check_complete().")

    def solve(self, assumptions=None):  # type: (Sequence[int]) -> bool
        """Solve the current set of clauses, optionally with a set of assumptions.

        Args:
            assumptions:
              An optional sequence of literals as integers, specified as in
              `add_clause()`.

        Returns:
            True if the clauses (and assumptions) are satisfiable, False otherwise.
        """
        if assumptions is None:
            return self.lib.solve(self.s)
        else:
            a = self._get_array(assumptions)
            a_ptr, size = self._to_intptr(a)
            return self.lib.solve_assumptions(self.s, size, a_ptr)

    def simplify(self):  # type: () -> bool
        '''Call Solver.simplify().'''
        return self.lib.simplify(self.s)

    def get_model(self, start=0, end=-1):  # type: (int, int) -> array.array
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

    def get_model_trues(self, start=0, end=-1, offset=0):  # type: (int, int, int) -> array.array
        """Get variables assigned true in the current model from the solver.

        Args:
            start, end (int):
              Optional start and end indices, interpreted as in ``range()``.
            offset (int):
              Optional offset to be added to the zero-based variable numbers
              from MiniSat.

        Returns:
            An array of true variables in the solver's current model.  If a
            start index was given, the variables are indexed from that value.
            """
        if end == -1:
            end = self.nvars()
        a = array.array('i', [-1] * (end-start))
        a_ptr, size = self._to_intptr(a)
        count = self.lib.getModelTrues(self.s, a_ptr, start, end, offset)
        # reduce the array down to just the valid indexes
        return a[:count]

    def block_model(self):
        """Block the current model from the solver."""
        model = self.get_model()
        self.add_clause([-(x+1) if model[x] > 0 else x+1 for x in range(len(model))])

    def model_value(self, i):  # type: (int) -> bool
        '''Get the value of a given variable in the current model.'''
        return self.lib.modelValue(self.s, i)

    def implies(self, assumptions=None):  # type(Sequence[int]) -> array.array
        """Get literals known to be implied by the current formula.  (I.e., all
        assignments made at level 0.)

        Args:
            assumptions:
              An optional sequence of literals as integers, specified as
              in `add_clause()`.

        Returns:
            An array of literals implied by the current formula (and optionally
            the given assumptions).
        """
        res = array.array('i', [-1] * self.nvars())
        res_ptr, _ = self._to_intptr(res)

        if assumptions is None:
            count = self.lib.getImplies(self.s, res_ptr)
        else:
            assumps = self._get_array(assumptions)
            assumps_ptr, assumps_size = self._to_intptr(assumps)
            count = self.lib.getImplies_assumptions(self.s, res_ptr, assumps_ptr, assumps_size)

        # reduce the array down to just the valid indexes
        return res[:count]

    def get_stats(self):
        """Returns a dictionary of solver statistics."""
        return {
            "solves": self.lib.get_solves(self.s),
            "starts": self.lib.get_starts(self.s),
            "decisions": self.lib.get_decisions(self.s),
            "rnd_decisions": self.lib.get_rnd_decisions(self.s),
            "propagations": self.lib.get_propagations(self.s),
            "conflicts": self.lib.get_conflicts(self.s),
        }


class SubsetMixin(Solver):
    """A mixin for any Solver class that lets it reason about subsets of a clause set."""
    _origvars = None    # type: int
    _relvars = None     # type: int

    def set_varcounts(self, vars, constraints):  # type: (int, int) -> None
        """Record how many of the solver's variables and clauses are
        "original," as opposed to clause-selector variables, etc.
        """
        self._origvars = vars
        self._relvars = constraints

    def add_clause_instrumented(self, lits, index):  # type: (Sequence[int], int) -> None
        """Add a "soft" clause with a relaxation variable (the relaxation var.
        is based on the index, which is assumed to be 0-based).

        Args:
            lits:
                A sequence of literals specified as in `add_clause()`.
            index (int):
                A 0-based index into the set of soft constraints.  The clause
                will be given a relaxation variable based on this index, and it
                will be used to specify the clause in subsets for
                `solve_subset()`, etc.
        """
        if self._origvars is None:
            raise Exception("SubsetSolver.set_varcounts() must be called before .add_clause_instrumented()")
        instrumented_clause = [-(self._origvars+1+index)]
        instrumented_clause.extend(lits)
        self.add_clause(instrumented_clause)

    def solve_subset(self, subset, extra_assumps=None):  # type: (Sequence[int], Sequence[int]) -> bool
        """Solve a subset of the constraints containing all "hard" clauses
        (those added with the regular `add_clause()` method) and the
        specified subset of soft constraints.

        Args:
            subset:
                A sequence of the indexes of any soft constraints to be included.
            extra_assumps:
                An optional sequence of extra literals to use when solving.

        Returns:
            True if the given subset is satisfiable, False otherwise.
        """
        if self._origvars is None:
            raise Exception("SubsetSolver.set_varcounts() must be called before .solve_subset()")

        assumptions = array.array('i', [i+self._origvars+1 for i in subset])
        if extra_assumps:
            assumptions.extend(extra_assumps)
        a_ptr, size = self._to_intptr(assumptions)
        return self.lib.solve_assumptions(self.s, size, a_ptr)

    def unsat_core(self, offset=0):  # type: (int) -> array.array
        """Get an UNSAT core from the last check performed by
        `solve_subset()`.  Assumes the last such check was UNSAT.

        Args:
            offset (int):
              Optional offset to be added to the zero-based indexes from
              MiniSat.

        Returns:
            An array of constraint indexes comprising an UNSAT core.
        """
        if self._origvars is None:
            raise Exception("SubsetSolver.set_varcounts() must be called (and at least one instrumented constraint added) before .unsat_core()")
        conflict_size = self.lib.conflictSize(self.s)
        a = array.array('i', [-1] * conflict_size)
        a_ptr, size = self._to_intptr(a)
        self.lib.unsatCore(self.s, self._origvars, a_ptr, offset)
        return a

    def sat_subset(self, offset=0):  # type: (int) -> array.array
        """Get the set of clauses satisfied in the last check performed by
        `solve_subset()`.  Assumes the last such check was SAT.  This may
        contain additional soft constraints not in the subset that was given to
        `solve_subset()`, if they were also satisfied by the model found.

        Args:
            offset (int):
              Optional offset to be added to the zero-based indexes from
              MiniSat.

        Returns:
            An array of constraint indexes comprising a satisfiable subset.
        """
        return self.get_model_trues(start=self._origvars, end=self._origvars+self._relvars, offset=offset)


class MinisatSolver(Solver):
    """A Python analog to MiniSat's Solver class.

    >>> S = MinisatSolver()

    Create variables using `new_var()`.  Add clauses as sequences of literals
    with `add_clause()`, analogous to MiniSat's ``addClause()``.  Literals are
    specified as integers, with the magnitude indicating the variable index
    (with 1-based counting) and the sign indicating True/False.  For example,
    to add clauses (x0), (!x1), (!x0 + x1 + !x2), and (x2 + x3):

    >>> for i in range(4):
    ...     S.new_var()  # doctest: +ELLIPSIS
    0
    1
    2
    3
    >>> for clause in [1], [-2], [-1, 2, -3], [3, 4]:
    ...     S.add_clause(clause)  # doctest: +ELLIPSIS
    True
    True
    True
    True

    The `solve()` method returns True or False just like MiniSat's.

    >>> S.solve()
    True

    Models are returned as arrays of Booleans, indexed by var.
    So the following represents x0=True, x1=False, x2=False, x3=True.

    >>> list(S.get_model())
    [1, 0, 0, 1]

    The `add_clause()` method may return False if a conflict is detected
    when adding the clause, even without search.

    >>> S.add_clause([-4])
    False
    >>> S.solve()
    False
    """
    def __init__(self):  # type: () -> None
        super(MinisatSolver, self).__init__("libminisat.so")


class MinicardSolver(Solver):
    """A Python analog to MiniCard's Solver class.

    >>> S = MinicardSolver()

    This has the same interface as `MinisatSolver`, with the addition of
    the `add_atmost()` and `add_atleast()` methods.

    >>> for i in range(4):
    ...     S.new_var()  # doctest: +ELLIPSIS
    0
    1
    2
    3
    >>> for clause in [1], [-2], [3, 4]:
    ...    S.add_clause(clause)
    True
    True
    True

    To add an AtMost constraint, specify the set of literals and the bound.  For example, to add AtMost({x0, !x1, x2}, 2):

    >>> S.add_atmost([1,-2,3], 2)
    True

    >>> S.solve()
    True

    >>> list(S.get_model())
    [1, 0, 0, 1]

    >>> S.add_atleast([1,2,3,4], 2)
    True

    >>> S.solve()
    True

    As with `add_clause()`, the `add_atmost()` method may return False if a
    conflict is detected when adding the constraint, even without search.

    >>> S.add_atmost([1,-3,4], 2)
    False
    >>> S.solve()
    False
    """
    def __init__(self):  # type: () -> None
        super(MinicardSolver, self).__init__("libminicard.so")

    def _setup_lib(self, libfilename):  # type: (str) -> None
        """Correct return types (if not int as assumed by ctypes) and set argtypes for
           functions from the minicard library.
        """
        super(MinicardSolver, self)._setup_lib(libfilename)

        # additional function for minicard
        l = self.lib
        l.addAtMost.restype = c_bool
        l.addAtMost.argtypes = [c_void_p, c_int, c_void_p, c_int]

    def add_atmost(self, lits, k):  # type: (Sequence[int], int) -> bool
        """Add an AtMost constraint to the solver.

        Args:
            lits:
              A sequence of literals as integers.  Each integer specifies a
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
            a = self._get_array(lits)
            a_ptr, size = self._to_intptr(a)
            return self.lib.addAtMost(self.s, size, a_ptr, k)
        else:
            return self.lib.addAtMost(self.s, 0, None, 0)

    def add_atleast(self, lits, k):  # type: (Sequence[int], int) -> bool
        """Convenience function to add an AtLeast constraint.
        Translates the AtLeast into an equivalent AtMost.
        See add_atmost().

        Args:
            lits:
              A sequence of literals as integers.  Each integer specifies a
              variable with **1**-based counting and a sign via the sign of
              the integer.  Ex.: [-1, 2, -3] is {!x0, x1, !x2}
            k (int):
              The [lower] bound to place on these literals.

        Returns:
            A boolean value returned from MiniCard's ``addAtMost()``
            function, indicating success (True) or conflict (False).
        """
        new_k = len(lits) - k
        new_lits = [-x for x in lits]
        return self.add_atmost(new_lits, new_k)


class MinisatSubsetSolver(SubsetMixin, MinisatSolver):
    """A class for reasoning about subsets of constraints within MiniSat.

    >>> S = MinisatSubsetSolver()

    It must be told explicitly how many of its variables are "real" and how
    many are relaxation variables for constraints.

    >>> S.set_varcounts(vars = 4, constraints = 5)

    And variables must be created for both the "real" variables and the
    relaxation variables.

    >>> for i in range(4+5):
    ...     _ = S.new_var()

    "Soft" clauses are added with `add_clause_instrumented()`, which has no
    return value, as it is impossible for these clauses to produce a conflict.

    >>> for i, clause in enumerate([[1], [-2], [-1, 2, 3], [-3], [-1]]):
    ...     S.add_clause_instrumented(clause, i)

    Any subset of the constraints can be tested for satisfiability.  Subsets
    are specified as sequences of soft clause indexes.

    >>> S.solve_subset([0,1,2])
    True

    Extra assumptions can be passed to solve_subset():

    >>> S.solve_subset([0,1,2], extra_assumps=[-3])
    False
    >>> S.solve_subset([0,1,2], extra_assumps=[3])
    True

    If a subset is found to be satisfiable, a potentially larger satisfied
    subset can be found.  Satisfiable subsets are returned as array objects.

    >>> satset = S.sat_subset()
    >>> sorted(satset)
    [0, 1, 2]

    If a subset is found to be unsatisfiable, an UNSAT core can be found.
    Cores are returned as array objects.

    >>> S.solve_subset([0,1,2,3])
    False

    >>> core = S.unsat_core()
    >>> sorted(core)
    [0, 1, 2, 3]
    """

    pass


class MinicardSubsetSolver(SubsetMixin, MinicardSolver):
    """A class for reasoning about subsets of constraints within MiniCard.

    This has the same interface as `MinisatSubsetSolver`, with the
    addition of the `add_atmost()` method.

    >>> S = MinicardSubsetSolver()
    >>> S.set_varcounts(vars = 4, constraints = 5)
    >>> for i in range(4+5):
    ...     _  = S.new_var()
    >>> for i, clause in enumerate([[1], [-2], [3], [4]]):
    ...     S.add_clause_instrumented(clause, i)

    AtMost and AtLeast constraints can be instrumented as well as clauses.

    >>> S.add_atmost_instrumented([1,-2,3], 2, 4)

    >>> S.solve_subset([0,1,4])
    True
    >>> S.solve_subset([0,1,2,3,4])
    False

    >>> core = S.unsat_core()
    >>> sorted(core)
    [0, 1, 2, 4]
    """

    def add_atmost_instrumented(self, lits, k, index):  # type: (Sequence[int], int, int) -> None
        """Add a "soft" ATMost constraint with a relaxation variable (the
        relaxation variable is based on the index, which is assumed to be
        0-based).

        Args:
            lits:
                A sequence of literals specified as in `add_atmost()`.
            k (int):
                The [upper] bound to place on these literals.
            index (int):
                A 0-based index into the set of soft constraints.  The clause
                will be given a relaxation variable based on this index, and it
                will be used to specify the clause in subsets for
                `solve_subset()`, etc.
        """
        if self._origvars is None:
            raise Exception("SubsetSolver.set_varcounts() must be called before .add_atmost_instrumented()")
        if not all(abs(x) <= self.nvars() for x in lits):
            raise Exception("Not all variables in %s are created yet.  Call new_var() first." % lits)
        if self._origvars+1+index > self.nvars():
            raise Exception("Relaxation variable %i has not been created yet.  Call new_var() first." % (self._origvars+1+index))

        numlits = len(lits)
        numnew = numlits - k
        # Original:     AtMost([lits], k)
        # Instrumented: AtMost([lits, newvar*numnew], numlits)
        # If newvar set False -> AtMost([lits], numlits) -> satisfied
        # If newvar set True  -> AtMost([lits], numlits-numnew) = AtMost([lits], k) [original]

        instrumented_lits = [self._origvars+1+index] * numnew
        instrumented_lits.extend(lits)
        a = self._get_array(instrumented_lits)
        a_ptr, size = self._to_intptr(a)
        self.lib.addAtMost(self.s, size, a_ptr, numlits)

    def add_atleast_instrumented(self, lits, k, index):  # type: (Sequence[int], int, int) -> None
        """Convenience function to add a soft AtLeast constraint.
        Translates the AtLeast into an equivalent AtMost.
        See add_atmost_instrumented().

        Args:
            lits:
                A sequence of literals specified as in `add_atmost()`.
            k (int):
                The [lower] bound to place on these literals.
            index (int):
                A 0-based index into the set of soft constraints.  The clause
                will be given a relaxation variable based on this index, and it
                will be used to specify the clause in subsets for
                `solve_subset()`, etc.
        """
        new_k = len(lits) - k
        new_lits = [-x for x in lits]
        return self.add_atmost_instrumented(new_lits, new_k, index)
