import numpy
import heapq
import cplex

from cplex.exceptions import CplexError

from Base import BiobjectiveSolver


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
            heapq.heappush(pq, (-RectangleSplittingSolver._calculate_rectangle_area(rectangle), rectangle, r_t, r_b))
            self._solutions.append(z_t)
            self._solutions.append(z_b)

            while pq:

                area, b, warm_t, warm_b = heapq.heappop(pq)
                print "\n\nCurrent Rectangle: ",b,"\n\n"
                rec_b = 0.5*(b[0][1]+b[1][1])

                #print "Current Box: ", b

                #search for new points in lower split
                z1_hat, r1_hat = self._lexmin(0, 1, rec_b, warmstart=warm_b)
                if not numpy.allclose(z1_hat.objs, b[1]):

                    #if found point is the same as initial point spanning the rectangle
                    #no unsupported point is in the rectangle -> discard the search
                    #else generate a new rectangle with the found point
                    #print "r1_hat ", z1_hat
                    rec1_hat = (z1_hat.objs, b[1])
                    self._solutions.append(z1_hat)
                    heapq.heappush(pq, (-RectangleSplittingSolver._calculate_rectangle_area(rec1_hat), rec1_hat, r1_hat, warm_b))

                #split the rectangle in vertical direction based on the new found point
                #and search in the upper half for new pareto points
                rec_t = z1_hat.objs[0]-BiobjectiveSolver.EPS
                z2_hat, r2_hat = self._lexmin(1, 0, rec_t, warmstart=warm_t,
                                              effort_level=0)
                if not numpy.allclose(z2_hat.objs, b[0]):

                    #again if the found point is the already known edge point
                    #one can discard the search in that rectangle
                    #print "r2_hat ", z2_hat
                    rec2_hat = (b[0], z2_hat.objs)
                    self._solutions.append(z2_hat)
                    heapq.heappush(pq, (-RectangleSplittingSolver._calculate_rectangle_area(rec2_hat),
                                        rec2_hat, warm_t, r2_hat))
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
                                init_rectangle, r_t, r_b))
            heapq.heappush(self._solutions, z_t)
            heapq.heappush(self._solutions, z_b)
            while pq:
                proof_empty = False
                area, b, warm_t, warm_b = heapq.heappop(pq)
                #print "AREA of current box", area
                #print "Current Rectangle: ",b
                rec_b = 0.5*(b[0][1]+b[1][1])

                #print "Current Box: ", b

                #search for new points in lower split
                z1_hat, r1_hat = self._lexmin(0, 1, rec_b, warmstart=warm_b)
                if not numpy.allclose(z1_hat.objs, b[1]):

                    #if found point is the same as initial point spanning the rectangle
                    #no unsupported point is in the rectangle -> discard the search
                    #else generate a new rectangle with the found point
                    #print "r1_hat ", z1_hat
                    rec1_hat = (z1_hat.objs, b[1])
                    heapq.heappush(self._solutions, z1_hat)
                    heapq.heappush(pq, (-RectangleSplittingSolver._calculate_rectangle_area(rec1_hat),
                                        rec1_hat, r1_hat, warm_b))
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
                z2_hat, r2_hat = self._lexmin(1, 0, rec_t, warmstart=warm_t, effort_level=0)
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
                    heapq.heappush(pq, (-RectangleSplittingSolver._calculate_rectangle_area(rec2_hat),
                                        rec2_hat, warm_t, r2_hat))
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
