import sys
import time

import cplex
from concurrent.rectangle_manager import RectangleSplittingManager
from algorithms.BiDirectionalEpsilonConstraint import DoubleEpsilonConstraintSolver
from algorithms.EpsilonConstraint import EpsilonConstraintSolver
from algorithms.RectangleSplitting import RectangleSplittingSolver
from algorithms.NormalConstraint import NormalConstraint

def rectangle():
    lp1 = sys.argv[1]
    lp2 = sys.argv[2]

    z1 = cplex.Cplex(lp1)
    z2 = cplex.Cplex(lp2)
    rectangle = RectangleSplittingSolver(z1, z2,  ["x", "y"])
    return rectangle.solve()


def epsilon():
    lp1 = sys.argv[1]
    lp2 = sys.argv[2]

    z1 = cplex.Cplex(lp1)
    z2 = cplex.Cplex(lp2)
    rectangle = EpsilonConstraintSolver(z1, z2,  ["x", "y"])
    return rectangle.solve()


def normal_constraint():
    lp1 = sys.argv[1]
    lp2 = sys.argv[2]

    z1 = cplex.Cplex(lp1)
    z2 = cplex.Cplex(lp2)
    nc = NormalConstraint(z1, z2, ["x", "y"])
    return nc.solve(4)


def double_epsilon():
    lp1 = sys.argv[1]
    lp2 = sys.argv[2]

    z1 = cplex.Cplex(lp1)
    z2 = cplex.Cplex(lp2)
    rectangle = DoubleEpsilonConstraintSolver(z1, z2,  ["x", "y"])
    return rectangle.solve()

    
def concurrent_rectangle():
    m = RectangleSplittingManager(sys.argv[1], sys.argv[2], ["x", "y"], 4)
    return m


if __name__ == "__main__":


    t = time.time()
    m = concurrent_rectangle()
    e = time.time()
    print "Construction of model Time: ", e-t

    t = time.time()
    sols = m.approximate(0.05)
    sols = normal_constraint()
    e = time.time()
    print "Concurrent New Rectangle Runtime: ", e-t
    for s in sols:
    #    #print s.vars
        print s.objs

    sys.exit()

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
