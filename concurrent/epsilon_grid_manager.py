from __future__ import division
import cplex
import numpy
import itertools

import multiprocessing as mp

from cplex.exceptions import CplexError

from concurrent.epsilon_grid_worker import EpsilonGridWorker
from concurrent.rectangle_worker import RectangleSplittingWorker
from utility.Hypervolume import HyperVolume


class EpsilonGridManager(object):
    """
        This is prarallel implementation of
        the epsilon contraint method with precomputed boundaries
    """

    def __init__(self, z1_name, z2_name, inter_vars, nof_worker):
        """
            init function

            @param z1_name: The ILP of the first objective
            @param z2_name: The ILP of the second objective
            @param nof_cpu: THe number of Processes to start
        """

        self.solutions = []
        self.models = (z1_name, z2_name)
        self.inter_vars = inter_vars
        self.nof_worker = nof_worker
        self.worker = []
        self._utopian_worker = []
        self._hypervol = HyperVolume()

        #concurrent stuff
        self.task_q = mp.JoinableQueue()
        self.task_utopian_q = mp.JoinableQueue()
        self.done_q = mp.Queue()

        z1_temp = cplex.Cplex(z1_name)
        z2_temp = cplex.Cplex(z2_name)
        #generate two worker dedicated to deal with the utopian point calculation
        for _ in xrange(2):
            z1 = cplex.Cplex(z1_temp)
            z2 = cplex.Cplex(z2_temp)
            #z1.parameters.mip.strategy.search.set(1)
            #z2.parameters.mip.strategy.search.set(1)
            z1.parameters.threads.set(max(int(mp.cpu_count()/2), 1))
            z2.parameters.threads.set(max(int(mp.cpu_count()/2), 1))
            z1.set_results_stream(None)
            z2.set_results_stream(None)

            p = RectangleSplittingWorker(z1, z2, ["z2_cons", "z1_cons"], inter_vars, self.task_utopian_q, self.done_q)
            p.deamon = True
            self._utopian_worker.append(p)
            p.start()

        #setup worker
        for _ in xrange(nof_worker):
            z1 = cplex.Cplex(z1_temp)
            z2 = cplex.Cplex(z2_temp)
            #z1.parameters.mip.strategy.search.set(1)
            #z2.parameters.mip.strategy.search.set(1)
            z1.parameters.threads.set(max(int(mp.cpu_count()/nof_worker), 1))
            z2.parameters.threads.set(max(int(mp.cpu_count()/nof_worker), 1))
            z1.set_results_stream(None)
            z2.set_results_stream(None)

            p = EpsilonGridWorker(z1, z2, inter_vars, self.task_q, self.done_q)
            p.deamon = True
            self.worker.append(p)
            p.start()

    def solve(self, nof_sol=None):

        nof_sol = len(self.worker) if nof_sol is None else nof_sol
        #init problems to solve
        self.task_utopian_q.put((0, 1, cplex.infinity, [None, None], ((None, None), (None, None))))
        self.task_utopian_q.put((1, 0, cplex.infinity, [None, None], ((None, None), (None, None))))

        self.task_utopian_q.join()

        b = [None, None]
        warm = [None, None]
        while not self.done_q.empty():
            pos, sol, warmstart, origin_rect = self.done_q.get()

            self.solutions.append(sol)
            b[pos] = sol.objs
            warm[pos] = warmstart[pos]

        delta = 1/(nof_sol - 1)
        diff = b[1][1] - b[0][1]
        alphas = [delta*i*diff for i in xrange(1, nof_sol - 1)]

        for i in xrange(nof_sol - 2):
            print "Bounds: ", b[0][1]+alphas[i]
            self.task_q.put(b[0][1]+alphas[i])

        self.task_q.join()

        while not self.done_q.empty():
            sols = self.done_q.get()
            self.solutions.append(sols)

        for _ in xrange(self.nof_worker):
            self.task_q.put("DONE")

        self.task_q.close()
        self.task_q.join()

        for p in self.worker:
            print "shut down worker ", p.pid
            p.join()

        for _ in xrange(2):
            self.task_utopian_q.put(["DONE", None, None, None, None])

        self.task_utopian_q.close()
        self.task_utopian_q.join()

        for p in self._utopian_worker:
            print "shut down worker ", p.pid
            p.join()

        return self.solutions

    def approximate(self, gap, nof_sol=None):
        return self.solve(nof_sol=nof_sol)