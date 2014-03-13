import sys
import time

import cplex
#from concurrent.concurrent_multiprocessing_cplex_old import *
from concurrent.rectangle_manager import RectangleSplittingManager
from single.multiobjective_cplex import *


#lp1 = sys.argv[1]
#lp2 = sys.argv[2]

#z1 = cplex.Cplex(lp1)
#z2 = cplex.Cplex(lp2)

def rectangle():
    lp1 = sys.argv[1]
    lp2 = sys.argv[2]

    z1 = cplex.Cplex(lp1)
    z2 = cplex.Cplex(lp2)
    rectangle = RectangleSplittingSolver(z1, z2, ["energy_cons", "immu_cons"], ["x", "y"])
    return rectangle.solve()

def epsilon():
    lp1 = sys.argv[1]
    lp2 = sys.argv[2]

    z1 = cplex.Cplex(lp1)
    z2 = cplex.Cplex(lp2)
    rectangle = EpsilonConstraintSolver(z1, z2, ["energy_cons", "immu_cons"], ["x", "y"])
    return rectangle.solve()

def double_epsilon():
    lp1 = sys.argv[1]
    lp2 = sys.argv[2]

    z1 = cplex.Cplex(lp1)
    z2 = cplex.Cplex(lp2)
    rectangle = DoubleEpsilonConstraintSolver(z1, z2, ["energy_cons", "immu_cons"], ["x", "y"])
    return rectangle.solve()

#def concurrent_rectangle_direct():
#    m =RectangleSplittingManager_directSplit(sys.argv[1], sys.argv[2], ["energy_cons", "immu_cons"], ["x", "y"], 4)
#    return m.solve()
    
def concurrent_rectangle():
    m = RectangleSplittingManager(sys.argv[1], sys.argv[2], ["energy_cons", "immu_cons"], ["x", "y"], 4)
    return m.solve()


#def concurrent_rectanlge_old():
#    m = RectangleSplittingManagerOld(sys.argv[1], sys.argv[2], ["energy_cons", "immu_cons"], ["x", "y"], 4)
#    return m.solve()


if __name__ == "__main__":
    #t=timeit.Timer("double_epsilon()", setup="from __main__ import double_epsilon")
    #print t.timeit(number=3)

    #cProfile.run("concurrent_rectangle()")

    #sys.exit()

   

    #t = time.time()
    #sols = concurrent_rectangle_direct()
    #e = time.time()
    #print "Concurrent New direct split Rectangle Runtime: ", e-t
    #for s in sols:
    #    #print s.vars
    #    print s.objs


    t = time.time()
    sols = concurrent_rectangle()
    e = time.time()
    print "Concurrent New Rectangle Runtime: ", e-t
    for s in sols:
    #    #print s.vars
        print s.objs

    sys.exit()

    t = time.time()
    sols = rectangle()
    e = time.time()
    print "Rectangle Runtime: ", e-t
    for s in sols:
    #    #print s.vars
        print s.objs

    t = time.time()
    sols = epsilon()
    e = time.time()
    print "Epsilon Runtime: ", e-t
    for s in sols:
    #    #print s.vars
        print s.objs



    t = time.time()
    sols = double_epsilon()
    e = time.time()
    print "Bidirectional Epsilon Runtime: ", e-t
    for s in sols:
    #    #print s.vars
        print s.objs


    sys.exit()


