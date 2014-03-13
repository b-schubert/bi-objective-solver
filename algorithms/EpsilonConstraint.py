import numpy
import heapq
import cplex

from cplex.exceptions import CplexError

from Base import BiobjectiveSolver


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
