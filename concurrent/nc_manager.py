from __future__ import division
import cplex
import numpy
import itertools

import multiprocessing as mp

from cplex.exceptions import CplexError

from concurrent.nc_worker import NormalConstraintWorker
from concurrent.rectangle_worker import RectangleSplittingWorker
from utility.Hypervolume import HyperVolume


class NormalConstraintManager(object):

    def __init__(self, z1_name, z2_name, inter_vars, nof_worker):
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
        m = mp.Manager()

        self.empty_rectangles = m.list()

        #generate two worker dedicated to deal with the utopian point calculation
        z_t = cplex.Cplex(z1_name)
        z_b = cplex.Cplex(z2_name)
        for _ in xrange(2):
            z1 = cplex.Cplex(z_t)
            z2 = cplex.Cplex(z_b)
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
            z1 = cplex.Cplex(z_t)
            z2 = cplex.Cplex(z_b)
            #z1.parameters.mip.strategy.search.set(1)
            #z2.parameters.mip.strategy.search.set(1)
            z1.parameters.threads.set(max(int(mp.cpu_count()/nof_worker), 1))
            z2.parameters.threads.set(max(int(mp.cpu_count()/nof_worker), 1))
            z1.set_results_stream(None)
            z2.set_results_stream(None)

            p = NormalConstraintWorker(z1, z2, inter_vars, self.task_q, self.done_q)
            p.deamon = True
            self.worker.append(p)
            p.start()

    def solve(self, nof_sol):

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

        print "Bound ",b
        z_t_bar = numpy.array(b[0])
        z_b_bar = numpy.array(b[1])
        print "Normalized obj, ", z_t_bar, z_b_bar
        #generate utopian line vector
        N1 = z_b_bar - z_t_bar
        print "N1 ", N1
        #generated normalized incremental step
        delta = 1/(nof_sol - 1)

        #generate utopian line points
        alphas = [delta*i for i in xrange(1, nof_sol - 1)]
        Xp = numpy.matrix([alphas[i]*z_t_bar + (1 - alphas[i])*z_b_bar for i in xrange(nof_sol - 2)])

        #generate tasks:
        for i in xrange(nof_sol - 2):
            self.task_q.put([Xp[i], N1, None])


        self.task_q.join()


        while not self.done_q.empty():
            sols = self.done_q.get()
            self.solutions.append(sols)



        for _ in xrange(self.nof_worker):
            self.task_q.put(["DONE", None, None])
        self.task_q.join()
        self.task_q.close()
        for p in self.worker:
            print "shut down worker ", p.pid
            p.join()

        for _ in xrange(2):
            self.task_utopian_q.put(["DONE", None, None, None, None])
        self.task_utopian_q.join()
        self.task_utopian_q.close()

        for p in self._utopian_worker:
            print "shut down worker ", p.pid
            p.join()

        return self.solutions