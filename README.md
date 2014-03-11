PyMiniSolvers
=============

PyMiniSolvers is a Python API to the MiniSat and MiniCard constraint solvers.
It closely matches the interfaces of their Solver classes, providing flexible
and powerful access to most of the solvers' standard capabilities.

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

