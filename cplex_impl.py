import sys
import time
import cPickle
import cplex
import re

from concurrent.rectangle_manager import RectangleSplittingManager
from concurrent.nc_manager import NormalConstraintManager
from concurrent.rectangle_nc_manager import RectangleNCManager

from algorithms.BiDirectionalEpsilonConstraint import DoubleEpsilonConstraintSolver
from algorithms.EpsilonConstraint import EpsilonConstraintSolver
from algorithms.RectangleSplitting import RectangleSplittingSolver
from algorithms.NormalConstraint import NormalConstraint

from utility.ParetoFilter import ParetoFilter
from utility.display import plot_pareto_front, plot_variable_distribution, plot_3d_variable_distribution

def hybrid():
    lp1 = sys.argv[1]
    lp2 = sys.argv[2]

    z1 = cplex.Cplex(lp1)
    z2 = cplex.Cplex(lp2)
    m = RectangleNCManager(z1, z2,  ["x", "y"], 4)
    return m


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
    nc = NormalConstraint(z1, z2, [])
    return nc.solve(15)


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

def concurrent_nc():
    m = NormalConstraintManager(sys.argv[1], sys.argv[2], ["x", "y"], 4)
    return m


if __name__ == "__main__":


    m = concurrent_rectangle()
    t = time.time()
    sols = sorted(m.approximate(0.0))
    e = time.time()
    print "Concurrent Rectangle Runtime: ", e-t
    for s in sols:
    #    #print s.vars
        print s.objs

    m = hybrid()
    t = time.time()
    sols = m.solve()
    e = time.time()
    print "Concurrent Hybrid Runtime: ", e-t
    for s in sols:
    #    #print s.vars
        print s.objs

    sys.exit()
    '''
    m = concurrent_nc()
    t = time.time()
    sols = m.solve(4)
    e = time.time()
    print "Concurrent NC Runtime: ", e-t
    for s in sols:
    #    #print s.vars
        print s.objs

    print
    filtered = ParetoFilter.filter(sols)
    for s in filtered:
    #    #print s.vars
        print s.objs

    sys.exit()

    t = time.time()
    sols = normal_constraint()
    e = time.time()
    print "New Normal Constraint Runtime: ", e-t
    for s in sols:
    #    #print s.vars
        print s.objs

    sys.exit()

    print
    filtered = ParetoFilter.filter(sols)
    for s in filtered:
    #    #print s.vars
        print s.objs

    sys.exit()

    m = concurrent_rectangle()
    t = time.time()
    sols = m.approximate(0.01)
    e = time.time()
    print "Concurrent New Rectangle Runtime: ", e-t

        t = time.time()
    sols = rectangle()
    e = time.time()
    print "Rectangle Runtime: ", e-t
    for s in sols:
    #    #print s.vars
        print s.objs
    sols.sort()
    '''




    m = concurrent_rectangle()
    t = time.time()
    sols = sorted(m.approximate(0.0))
    e = time.time()
    print "Concurrent New Rectangle Runtime: ", e-t

    plot_pareto_front("../test_pareto_front.png", map(lambda x: x.objs[0],sols),map(lambda x: x.objs[1],sols))


    c = 0
    xes = []
    yes = []
    for s in sols:
    #    #print s.vars
        print s.objs


        seq = filter(lambda x: "x" in x, s.vars)

        #print seq
        epi = filter(lambda x: "y" in x[0], s.vars.items())
        #print epi
        wt = "EVRNPCACFRNYVPVCGSDGKTYGNPCMLNCAAQTKVPGLKLVHEGRCQRSNVEQF"
        mutant = ["A"]*len(wt)
        ys = [0]*len(wt)
        epitopes = {}
        for x in seq:
            pos,aa =  re.match("x\((.*),(.*)\)",x).group(1,2)
            mutant[int(pos)-1] = aa
        for y in epi:
            print y, y[0],y[1]
            allele,pos =  re.match("y\((.*),(.*)\)",y[0]).group(1,2)
            pos = int(pos)-1
            epitopes.setdefault(allele,[]).append(pos)
            for i in xrange(9):
                ys[pos+i] += float(y[1])

        print "Ys", ys
        xes.append(list(mutant))
        yes.append(ys)
        seq_len = len(wt)
        plot_variable_distribution("../test_variables_correlation_%i.png"%c, list(mutant), ys)
        print "mt       : "+"".join(mutant[i] for i in xrange(seq_len))
        c+=1
        seq_len = len(wt)
        print "Hemming Distance: "+str(sum( wt[i] != mutant[i] for i in range(seq_len) ))+ str([ (i+1,wt[i],mutant[i]) for i in range(seq_len) if wt[i] != mutant[i]])
        print "wt       : "+wt
        print "mt       : "+"".join(mutant[i] for i in xrange(seq_len))

        for allele, ep_pos in sorted(epitopes.items()):
                    st = allele+": "
                    i=0

                    while i < seq_len:
                        if i in epitopes[allele]:
                            end = i+9 if i+9 <= seq_len else seq_len
                            st+="".join( "^" if wt[j] != mutant[j] else "*"  for j in xrange(i,end) )
                            i+=9
                        elif wt[i] != mutant[i]:
                            st+="^"
                            i+=1
                        else:
                            st+="-"
                            i+=1
                    print st
    plot_3d_variable_distribution("../test_3d_var_dist.png", xes, yes)
    sys.exit()


    m = concurrent_rectangle()
    t = time.time()
    sols = m.approximate(0.01)
    e = time.time()
    print "Concurrent New Rectangle Runtime: ", e-t



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
