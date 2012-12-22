import minisat
import minicard

for s in [minisat.Solver(), minicard.Solver()]:

    for i in range(5):
        s.new_var()
        print s.nvars()

    s.add_clause([-1,2])
    s.add_clause([1,-2])
    s.add_clause([-1,-2,3,4,5])
    s.add_clause([1])

    print "---"
    print s.nclauses()
    s.simplify()
    print s.nclauses()
    sat = s.solve()
    print "---"
    if sat:
        print "sat"
        for i in range(5):
            print s.model_value(i)
    else:
        print "unsat"

