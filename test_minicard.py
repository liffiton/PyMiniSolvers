
import minicard
import unittest

class MinicardTest(unittest.TestCase):
    def setUp(self):
        self.solver = minicard.Solver()
        self.atmosts = [(list(range(1,31)), 15)]
        self.clauses = [ [15], [10], [8], [12, 5], [4, 6], [2], [1], [-1, -2, -6, -8, -10, -12, -15] ]
        self.assumptions = [-x for x in range(16,31)]
        self.numvars = 30

    def tearDown(self):
        del self.solver

    def test_newvars(self):
        for i in range(self.numvars):
            self.solver.new_var(True)
            self.assertEqual(self.solver.nvars(), i+1)

    def test_add_clause_without_vars(self):
        self.assertRaises(Exception, self.solver.add_clause, [-1,2])

    def make_vars(self):
        for i in range(self.numvars):
            self.solver.new_var()

    def add_subset(self, subset):
        for cl in subset:
            self.solver.add_clause(cl)

    def add_atmosts(self, atmosts):
        for atmost in atmosts:
            self.solver.add_atmost(atmost[0], atmost[1])

    def test_sat(self):
        self.make_vars()
        self.add_subset(self.clauses[:-1])
        self.assertEqual(self.solver.solve(), True)

    def test_unsat(self):
        self.make_vars()
        self.add_subset(self.clauses)
        self.solver.add_clause([-8, -10])
        self.assertEqual(self.solver.solve(), False)

    def test_assumptions(self):
        self.make_vars()
        self.add_subset(self.clauses[:-2])
        self.assertEqual(self.solver.solve([-5]), True)
        self.assertEqual(self.solver.solve([-6]), True)
        self.assertEqual(self.solver.solve([-8,-10]), False)

    def test_model(self):
        self.make_vars()
        from math import copysign
        isPositive = lambda x: copysign(1,x) > 0
        subset = self.clauses[:-1]
        self.add_subset(subset)
        self.solver.solve()
        m = self.solver.get_model()
        self.assertEqual(len(m), self.solver.nvars())
        for cl in subset:
            self.assertTrue(any([ m[abs(x)-1] == isPositive(x) for x in cl ]))

    def int_check(self):
        import random
        for i in range(1000):
            assumps = [x*random.choice([1,-1]) for x in self.assumptions]
            self.solver.solve(assumps)
            self.assertEqual(self.solver.solve(self.assumptions), True)

    def test_complete(self):
        self.make_vars()
        self.add_atmosts(self.atmosts)
        self.int_check()
        self.add_subset(self.clauses)
        self.int_check()
        self.assertEqual(self.solver.solve(self.assumptions), True)

if __name__ == '__main__':
    unittest.main()

