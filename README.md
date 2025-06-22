PyMiniSolvers
=============

PyMiniSolvers is a Python API for the [MiniSat](http://minisat.se/) and
[MiniCard](http://git.io/minicard) constraint solvers.  It accesses the solvers
via shared libraries, reducing overhead and allowing for efficient incremental
solving.  Its interfaces closely match the interfaces of each tool's Solver
class, providing flexible and powerful access to most of the solvers' standard
capabilities.  Additionally, an extended "SubsetSolver" interface is provided
for each tool, simplifying reasoning about subsets of a constraint set (e.g.,
solving arbitrary subsets of the added constraints, extracting UNSAT cores,
etc.).

Basic usage:
```python
>>> import minisolvers
>>> S = minisolvers.MinisatSolver()
>>> S.new_vars(4)  
>>> for clause in [1], [-2], [-1, 2, -3], [3, 4]:
...     S.add_clause(clause)  

>>> S.solve()
True

>>> list(S.get_model())
[1, 0, 0, 1]

>>> S.add_clause([-1, 2, 3, -4])

>>> S.solve()
False
```

For further examples, see the documentation.

API Documentation
-----------------

https://pyminisolvers.readthedocs.io/

Setup
-----

Requirements:
 - Python 3.9+
 - A standard build environment (make, gcc, etc.)

To build the shared libraries:

    $ make

To test the API:

    $ make test

Or manually run the tests:

    $ python3 -m doctest -v minisolvers.py
    $ python3 test_minisolvers.py

License
-------

This code is licensed under the MIT license.  See LICENSE for details.

