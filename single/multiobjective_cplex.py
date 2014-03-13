"""
This file contains biobjective solver written in
CPLEX Python API

@Author: Benjamin Schubert
@Version: 0.1
"""
import numpy
import heapq
import itertools
import cplex
import abc

from cplex.exceptions import CplexError

from utility.Solution import Solution
from utility.Hypervolume import HyperVolume


class BiobjectiveSolver(object):
    """
        Generic biobjectiv solver
    """

    __metaclass__ = abc.ABCMeta

    EPS = 1e-2
    #EPS = 1

    def __init__(self, z1, z2, biob_constraints, interesting_vars):
        """
            initializes the solver model and
             modifies the given cplex model

             @param z1: cplex model with first objective
             @param z2: equivalent cplex model with second objective
             @param biob_constraints: a list of constraints which have to be changed in each iteration
                    should be only two!
             @param interesting_vars: list of variable names of the model which a user is interested in
        """
        assert len(biob_constraints) == 2
        #debug
        z1.set_results_stream(None)
        z2.set_results_stream(None)

        self._models = (z1, z2)
        self._changeable_constraints = ["z2_cons", "z1_cons"]
        self._inter_variables = filter(lambda x:  x[0] in interesting_vars, z1.variables.get_names())
        self._solutions = []
        self._variables = z1.variables.get_names()
        self._hypervol = HyperVolume()

        #transform problems by adding an additional box constraint to each of the two
        #concerning the other objective function

        #dont know why i first have to go through the complete list and filter it
        z1_obj_val = { v:z1.objective.get_linear(v) for v in self._variables if z1.objective.get_linear(v) != 0.0}
        z2_obj_val = { v:z2.objective.get_linear(v) for v in self._variables if z2.objective.get_linear(v) != 0.0}
        z1.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=z2_obj_val.keys(), val=z2_obj_val.values())], senses=["L"], rhs=[0.0], range_values=[0], names=["z2_cons"])
        z2.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=z1_obj_val.keys(), val=z1_obj_val.values())], senses=["L"], rhs=[0.0], names=["z1_cons"])

    def _lexmin(self, z1_idx, z2_idx, boundary,  warmstart=None, effort_level=0):
        """
            lexicographic optimization based on input values

            @param z1_idx, z2_idx: integers defining which of the two objectives is optimized first
            @param z1_bound, boundary: floating point of new boundary of the
            @param warmstart (optional): is an optional feature to use a warmstart for z1
            @return Solution, values_of_all_variables_for_warm_start


            @QUESTION:Dont know if one should first delete old MIP starts or just keep adding them
        """

        z1 = self._models[z1_idx]
        if warmstart:
            z1.MIP_starts.add([self._variables, warmstart], effort_level)
        z1.linear_constraints.set_rhs(self._changeable_constraints[z1_idx], boundary)
        z1.solve()
        z1_hat = z1.solution.get_objective_value()
        z1_hat_values = z1.solution.get_values()

        #lexicographical solution second model
        z2 = self._models[z2_idx]

        z2.linear_constraints.set_rhs(self._changeable_constraints[z2_idx], z1_hat)
        #set warm start from first objective
        z2.MIP_starts.add([self._variables, z1_hat_values], z2.MIP_starts.effort_level.auto)
        z2.solve()
        z2_hat = z2.solution.get_objective_value()
        z2_hat_values = z2.solution.get_values()

        #generate solution object
        objs = [0, 0]
        objs[z1_idx] = z1_hat
        objs[z2_idx] = z2_hat

        inter_vars = {}
        if len(self._inter_variables) > 0:
            inter_vars = {k:v for v, k in itertools.izip(z2.solution.get_values(self._inter_variables), self._inter_variables)
                          if v > 0.0}
        s = Solution(objs, inter_vars)

        return s, z2_hat_values

    @abc.abstractmethod
    def solve(self):
        """
            solver functionality
        """
        NotImplementedError("solve function has to be implemented")

    @abc.abstractmethod
    def approximate(self, gap):
        """
            implements an approximation technique
        """
        NotImplementedError("solve function has to be implemented")

    def get_solutions(self):
        return self._solutions


class RectangleSplittingSolver(BiobjectiveSolver):
    """
        Rectangle splitting algorithm for bi-objective MIPs based on
        Noland et al., Criterion Space Search Algorithm for MIP,  Nov 2013, Optimization Online
    """


    @staticmethod
    def _calculate_rectangle_area(b):
        """
            Calculates the area of a
        """
        yd = b[0][1] - b[1][1]
        xd = b[0][0] - b[1][0]
        return numpy.sqrt(yd*yd)+numpy.sqrt(xd*xd)

    def solve(self):
        """
            implements the Rectangle Splitting Algorithm
        """
        #Priority Queue: prioritization based on area of the search rectangle
        self._solutions = []
        pq = []
        #generate initial vertices of the rectangle
        try:
            z_t, r_t = self._lexmin(0, 1, cplex.infinity)
            z_b, r_b = self._lexmin(1, 0, cplex.infinity)
            #spann the rectangle
            rectangle = (z_t.objs, z_b.objs)
            #print "Start Rectangle: ",rectangle
            heapq.heappush(pq, (-RectangleSplittingSolver._calculate_rectangle_area(rectangle), rectangle, r_b))
            self._solutions.append(z_t)
            self._solutions.append(z_b)

            while pq:

                area, b, warmstart = heapq.heappop(pq)
                print "\n\nCurrent Rectangle: ",b,"\n\n"
                rec_b = 0.5*(b[0][1]+b[1][1])

                #print "Current Box: ", b

                #search for new points in lower split
                z1_hat, r1_hat = self._lexmin(0, 1, rec_b, warmstart=warmstart)
                if not numpy.allclose(z1_hat.objs, b[1]):

                    #if found point is the same as initial point spanning the rectangle
                    #no unsupported point is in the rectangle -> discard the search
                    #else generate a new rectangle with the found point
                    #print "r1_hat ", z1_hat
                    rec1_hat = (z1_hat.objs, b[1])
                    self._solutions.append(z1_hat)
                    heapq.heappush(pq, (-RectangleSplittingSolver._calculate_rectangle_area(rec1_hat), rec1_hat, warmstart))

                #split the rectangle in vertical direction based on the new found point
                #and search in the upper half for new pareto points
                rec_t = z1_hat.objs[0]-BiobjectiveSolver.EPS
                z2_hat, r2_hat = self._lexmin(1, 0, rec_t, warmstart=r1_hat,
                                              effort_level=4)
                if not numpy.allclose(z2_hat.objs, b[0]):

                    #again if the found point is the already known edge point
                    #one can discard the search in that rectangle
                    #print "r2_hat ", z2_hat
                    rec2_hat = (b[0], z2_hat.objs)
                    self._solutions.append(z2_hat)
                    heapq.heappush(pq, (-RectangleSplittingSolver._calculate_rectangle_area(rec2_hat), rec2_hat, r2_hat))
        except CplexError, exc:
            print exc
            return

        return self._solutions

    def approximate(self, eps=0.05):
        """
            approximates the pareto front with an epsilon bound on precission.
            it uses the notion of hypervolume indicator (resp. modified hypervolume indicator)

            it calculates the area of the intersecting polygon of the currently found pareto points,
            (zt,zb)- R^2_>= [ H(y)], and the spanning search rectangles [h(y)]. The differ h(y) - H(y)
            is the quantitative measure of the approximation of the pareto front

            for calculating the areas a clipping algorithm has been used (don't know jet which one)

            @param eps: the precision of approximation [0.0,1.0] where 0.0 indicates an exact calculation
            @return: A list of pareto points as Solution objects
        """

        if eps <= 0.0:
            return self.solve()

        #here starts the actual work. One has to generate different structures. you have to know:
        #1) The initial search rectangle (this is the area on which we will clip)
        #2) All search rectangles that the current pareto points span
        #3) If for a search rectangle one has proven that no further pareto point lies within it
        #   (those rectangles can be excluded from the list of search rectangles used for the hypervolume calculation)
        self._solutions = []
        pq = []

        #this will save if we could prove if in a certain rectangle is no other pareto point
        search_rectangles = set()

        try:
            z_t, r_t = self._lexmin(0, 1, cplex.infinity)
            z_b, r_b = self._lexmin(1, 0, cplex.infinity)
            #spann the rectangle
            init_rectangle = (z_t.objs, z_b.objs)
            #print "Start Rectangle: ",rectangle
            heapq.heappush(pq, (-RectangleSplittingSolver._calculate_rectangle_area(init_rectangle),
                                init_rectangle, r_b))
            heapq.heappush(self._solutions, z_t)
            heapq.heappush(self._solutions, z_b)
            while pq:
                proof_empty = False
                area, b, warmstart = heapq.heappop(pq)
                #print "AREA of current box", area
                #print "Current Rectangle: ",b
                rec_b = 0.5*(b[0][1]+b[1][1])

                #print "Current Box: ", b

                #search for new points in lower split
                z1_hat, r1_hat = self._lexmin(0, 1, rec_b, warmstart=warmstart)
                if not numpy.allclose(z1_hat.objs, b[1]):

                    #if found point is the same as initial point spanning the rectangle
                    #no unsupported point is in the rectangle -> discard the search
                    #else generate a new rectangle with the found point
                    #print "r1_hat ", z1_hat
                    rec1_hat = (z1_hat.objs, b[1])
                    heapq.heappush(self._solutions, z1_hat)
                    heapq.heappush(pq, (-RectangleSplittingSolver._calculate_rectangle_area(rec1_hat),
                                        rec1_hat, warmstart))
                else:
                    proof_empty = True

                #dont know if one should check every time....
                gap = self._hypervol.calc_hypervol_gap(self._solutions, init_rectangle, search_rectangles)
                #print "Approximation Gap: ", gap
                #print "Proven empty Rectangles ", search_rectangles
                if gap <= eps:
                    return self._solutions

                #split the rectangle in vertical direction based on the new found point
                #and search in the upper half for new pareto points
                rec_t = z1_hat.objs[0]-BiobjectiveSolver.EPS
                z2_hat, r2_hat = self._lexmin(1, 0, rec_t, warmstart=r1_hat, effort_level=4)
                if not numpy.allclose(z2_hat.objs, b[0]):

                    #again if the found point is the already known edge point
                    #one can discard the search in that rectangle, because
                    #here one has proven that there cannot be further points
                    #within the rectangle
                    #print "r2_hat ", z2_hat

                    #if z1_hat could not be found -> (z2_hat,z_b) is empty
                    #if z1_hat could be found -> (z2_hat,z1_hat) is empty
                    if proof_empty:
                        search_rectangles.add((z2_hat.objs, b[1]))
                    else:
                        search_rectangles.add((z2_hat.objs, z1_hat.objs))
                    rec2_hat = (b[0], z2_hat.objs)
                    heapq.heappush(self._solutions, z2_hat)
                    heapq.heappush(pq, (-RectangleSplittingSolver._calculate_rectangle_area(rec2_hat), rec2_hat, r2_hat))
                else:
                    #if z1_hat could not be found initial rectangle was empty -> b
                    #if z1_hat was found but z2_hat was not (z_t,z1_hat) is empty
                    if proof_empty:
                        search_rectangles.add(b)
                    else:
                        search_rectangles.add((b[0], z1_hat.objs))

                gap = self._hypervol.calc_hypervol_gap(self._solutions, init_rectangle, search_rectangles)
                #print "Approximation Gap: ", gap
                #print "Proofed empty Rectangles ", search_rectangles
                if gap <= eps:
                    return self._solutions

        except CplexError, exc:
            print exc
            return

        return self._solutions


class EpsilonConstraintSolver(BiobjectiveSolver):
    """
        Implements the naive Epsilon Constraint approach
        that is lexmin(z1,z2, R(z^l -(0,EPS)',z^B) until z^B is reached
    """

    def solve(self):
        """
            solves the epsilon constraint biobjective problem
        """
        pq = []
        #generate initial vertices of the rectangle
        try:
            z_t, r_t = self._lexmin(0, 1, cplex.infinity)
            z_b, r_b = self._lexmin(1, 0, cplex.infinity)
            print z_t,z_b
            #print "Lower Vertex ", z_b.objs
            self._solutions.append(z_t)
            self._solutions.append(z_b)
            pq.append(z_t.objs[1])

            while pq:
                bound = pq.pop() - BiobjectiveSolver.EPS
                z1_hat, r_hat = self._lexmin(0, 1, bound)
                print "\n\nCurrent bound", bound, z1_hat.objs,"\n\n"
                if not numpy.allclose(z1_hat.objs, z_b.objs):
                    #if we have not reached the bottom corner of the search space add
                    #add new boundary and continue
                    pq.append(z1_hat.objs[1])
                    self._solutions.append(z1_hat)
        except CplexError, exc:
            print exc
            return

        return self._solutions

    def approximate(self, eps):
        """
            approximates the pareto front with an epsilon bound on precission.
            it uses the notion of hypervolume indicator (resp. modified hypervolume indicator)

            it calculates the area of the intersecting polygon of the currently found pareto points,
            (zt,zb)- R^2_>= [ H(y)], and the spanning search rectangles [h(y)]. The differ h(y) - H(y)
            is the quantitative measure of the approximation of the pareto front

            for calculating the areas a clipping algorithm has been used (don't know jet which one)

            @param eps: the precision of approximation [0.0,1.0] where 0.0 indicates an exact calculation
            @return: A list of pareto points as Solution objects
        """
        if eps <= 0.0:
            return self.solve()

        pq = []
        #generate initial vertices of the rectangle
        try:
            z_t, r_t = self._lexmin(0, 1, cplex.infinity)
            z_b, r_b = self._lexmin(1, 0, cplex.infinity)

            #print "Lower Vertex ", z_b.objs
            init_rectangle = (z_t.objs, z_b.objs)
            heapq.heappush(self._solutions, z_t)
            heapq.heappush(self._solutions, z_b)
            pq.append(z_t.objs)

            search_rectangles = set()
            while pq:
                gap = self._hypervol.calc_hypervol_gap(self._solutions, init_rectangle, search_rectangles)
                #print "Hypergap ", gap
                if gap <= eps:
                    return self._solutions

                point = pq.pop()
                bound = point[1] - BiobjectiveSolver.EPS
                z1_hat, r_hat = self._lexmin(0, 1, bound)
                #print bound, z1_hat.objs
                if not numpy.allclose(z1_hat.objs, z_b.objs):
                    #if we have not reached the bottom corner of the search space add
                    #add new boundary and continue
                    rectangle = (point, z1_hat.objs)
                    self._solutions.append(z1_hat)
                    pq.append(z1_hat.objs)
                    search_rectangles.add(rectangle)
        except CplexError, exc:
            print exc
            return

        return self._solutions


class DoubleEpsilonConstraintSolver(BiobjectiveSolver):
    """
        Solves a modification of the epsilon constraint problem.
        Epsilon bound is applied from top and bottom vertices in
        an alternating fashion, such that the previous solved
        problem can be used as warm start to the current one
    """

    def solve(self):
        """
            solves the double epsilon constraint problem and returns the exact pareto front
        """
        pq = []
        try:
            #init problem/rectangle
            z_t, r_t = self._lexmin(0, 1, cplex.infinity)
            z_b, r_b = self._lexmin(1, 0, cplex.infinity)
            self._solutions.append(z_t)
            self._solutions.append(z_b)
            pq.append((z_t.objs, z_b.objs, r_b))

            while pq:
                z1_objs, z2_objs, warmstart = pq.pop()
                print "\n\nCurrent bounds" , z1_objs, z2_objs,"\n\n"
                z1_hat, r1_hat = self._lexmin(0, 1, z1_objs[1] - BiobjectiveSolver.EPS, warmstart=warmstart)

                if numpy.allclose(z1_hat.objs, z2_objs):
                    break

                z2_hat, r2_hat = self._lexmin(1, 0, z2_objs[0] - BiobjectiveSolver.EPS, warmstart=r1_hat)

                if numpy.allclose(z2_hat.objs, z1_objs) or numpy.allclose(z2_hat.objs, z1_hat.objs):
                    break

                self._solutions.append(z1_hat)
                self._solutions.append(z2_hat)
                pq.append((z1_hat.objs, z2_hat.objs, r2_hat))

        except CplexError, exc:
            print exc
            return

        return self._solutions

    def approximate(self, eps):
        """
            approximates the pareto front with an epsilon bound on precission.
            it uses the notion of hypervolume indicator (resp. modified hypervolume indicator)

            it calculates the area of the intersecting polygon of the currently found pareto points,
            (zt,zb)- R^2_>= [ H(y)], and the spanning search rectangles [h(y)]. The differ h(y) - H(y)
            is the quantitative measure of the approximation of the pareto front

            for calculating the areas a clipping algorithm has been used (don't know jet which one)

            @param eps: the precision of approximation [0.0,1.0] where 0.0 indicates an exact calculation
            @return: A list of pareto points as Solution objects
        """
        if eps <= 0.0:
            return self.solve()

        pq = []
        search_rectangles = set()
        try:
            #init problem/rectangle
            z_t, r_t = self._lexmin(0, 1, cplex.infinity)
            z_b, r_b = self._lexmin(1, 0, cplex.infinity)
            self._solutions.append(z_t)
            init_rectangle = (z_t.objs, z_b.objs)
            self._solutions.append(z_b)

            pq.append((z_t.objs, z_b.objs, r_b))

            while pq:
                z1_objs, z2_objs, warmstart = pq.pop()
                #print z1_objs, z2_objs,
                z1_hat, r1_hat = self._lexmin(0, 1, z1_objs[1] - BiobjectiveSolver.EPS, warmstart=warmstart)
                self._solutions.append(z1_hat)
                search_rectangles.add((z1_objs, z1_hat.objs))

                gap = self._hypervol.calc_hypervol_gap(self._solutions, init_rectangle, search_rectangles)
                #print "Hypergap ", gap
                if gap <= eps:
                    return self._solutions

                if numpy.allclose(z1_hat.objs, z2_objs):
                    return self._solutions

                z2_hat, r2_hat = self._lexmin(1, 0, z2_objs[0] - BiobjectiveSolver.EPS, warmstart=r1_hat)
                self._solutions.append(z2_hat)
                search_rectangles.add((z2_hat.objs, z2_objs))

                gap = self._hypervol.calc_hypervol_gap(self._solutions, init_rectangle, search_rectangles)
                #print "Hypergap ", gap
                if gap <= eps:
                    return self._solutions

                if numpy.allclose(z2_hat.objs, z1_objs) or z2_hat == z1_hat:
                    return self._solutions

                pq.append((z1_hat.objs, z2_hat.objs, r2_hat))

        except CplexError, exc:
            print exc
            return

        return self._solutions


class NormalConstraint(BiobjectiveSolver):
    """
        Implements the normal constraint method of
        Messac et al., 2003, "The Normalized Normal Constraint Method for Generating the Pareto Frontier",
        Structural and Multidisciplinary Optimization, Vol. 25, No. 2
        for a bi-objective case
    """

    def __init__(self, z1, z2, biob_constraints, interesting_vars):
        super(NormalConstraint, self).__init__(z1, z2, biob_constraints, interesting_vars)