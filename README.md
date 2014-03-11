PyMiniSolvers
=============

PyMiniSolvers is a Python API for the [MiniSat](http://minisat.se/) and
[MiniCard](http://git.io/minicard) constraint solvers.  It closely matches the
interfaces of their Solver classes, providing flexible and powerful access to
most of the solvers' standard capabilities.

Example usage:
```python
>>> import minisolvers
>>> S = minisolvers.MinisatSolver()
>>> for i in range(4):
...     S.new_var()  
>>> for clause in [1], [-2], [-1, 2, -3], [3, 4]:
...     S.add_clause(clause)  

>>> S.solve()
True

>>> list(S.get_model())
[1, 0, 0, 1]
```

For further examples, see the documentation.

API Documentation
-----------------

http://liffiton.github.io/PyMiniSolvers/

Setup
-----

Requirements:
 - Python 2.7 or above (compatible with Python 3)
 - A standard build environment (make, gcc, etc.)

Tested Platforms:
 - Linux
 - Cygwin
 - OS X

To build and test the API:

    $ make
    $ python -m doctest -v minisolvers.py
    $ python test_minisolvers.py

License
-------

This code is licensed under the MIT license.  See LICENSE for details.

