import minisat
import unittest

class MinisatTest(unittest.TestCase):
    def setUp(self):
        self.solver = minisat.Solver()
        self.clauses = [ [1], [-2], [3, 4], [-3, 5], [-4, 6], [-5, 4], [-6] ]
        self.numvars = max([max(cl) for cl in [[abs(x) for x in cl] for cl in self.clauses]])

    def tearDown(self):
        del self.solver

    def test_newvars(self):
        for i in range(self.numvars):
            self.solver.new_var()
            self.assertEqual(self.solver.nvars(), i+1)

    def test_add_clause_without_vars(self):
        self.assertRaises(Exception, self.solver.add_clause, [-1,2])

    def add_subset(self, subset):
        for i in range(self.numvars):
            self.solver.new_var()
        for cl in subset:
            self.solver.add_clause(cl)

    def test_sat(self):
        self.add_subset(self.clauses[:-1])
        self.assertEqual(self.solver.solve(), True)

    def test_unsat(self):
        self.add_subset(self.clauses)
        self.assertEqual(self.solver.solve(), False)

    def test_assumptions(self):
        self.add_subset(self.clauses[:-2])
        self.assertEqual(self.solver.solve([-5]), True)
        self.assertEqual(self.solver.solve([-6]), True)
        self.assertEqual(self.solver.solve([-5,-6]), False)

    def test_model(self):
        from math import copysign
        isPositive = lambda x: copysign(1,x) > 0
        subset = self.clauses[:-1]
        self.add_subset(subset)
        self.solver.solve()
        m = self.solver.get_model()
        self.assertEqual(len(m), self.solver.nvars())
        for cl in subset:
            self.assertTrue(any([ m[abs(x)-1] == isPositive(x) for x in cl ]))


if __name__ == '__main__':
    unittest.main()

