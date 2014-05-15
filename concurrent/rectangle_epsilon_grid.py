from __future__ import division
import cplex
import numpy
import itertools

import multiprocessing as mp

from cplex.exceptions import CplexError

from concurrent.epsilon_grid_manager import EpsilonGridManager
from concurrent.rectangle_manager import RectangleSplittingManager

from utility.Hypervolume import HyperVolume
from utility.ParetoFilter import ParetoFilter


class RectangleEpsilonGridManager(object):

    def __init__(self, z1_name, z2_name, inter_vars, nof_worker):
        self.solutions = []
        self._nc = EpsilonGridManager(z1_name, z2_name, inter_vars, nof_worker+2)
        print "Epsilon Grid initialized"
        self._rs = RectangleSplittingManager(z1_name, z2_name, inter_vars, nof_worker)
        print "Rectangle Splitting initialized"
        self._nof_worker = nof_worker

    def solve(self):

        sols = self._nc.solve(self._nof_worker)
        print "unfiltered pareto points"
        for s in sols:
            print s

        #filter non pareto points
        sols = sorted(ParetoFilter.filter(sols))
        #self.solutions.extend(sols)

        print "Filtered pareto points"
        for s in sols:
            print s

        sols = self._rs.solve(init_recs=sols)
        self.solutions.extend(sols)
        return self.solutions

    def approximate(self, gap):
        sols = self._nc.solve(self._nof_worker)

        #filter non pareto points
        sols = sorted(ParetoFilter.filter(sols))
        #self.solutions.extend(sols)

        print "Filtered pareto points"
        for s in self.solutions:
            print s

        curr_gap = HyperVolume.calc_hypervol_gap(self.solutions, [sols[0].objs, sols[-1].objs], [])
        if numpy.allclose(gap, curr_gap, rtol=1e-01, atol=1e-04) or curr_gap < gap:
            return self.solutions

        sols = self._rs.approximate(gap, init_recs=sols)
        self.solutions.extend(sols)

        return self.solutions