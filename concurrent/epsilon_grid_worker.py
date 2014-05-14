import numpy
import cplex
import abc
import itertools
import multiprocessing as mp

from cplex.exceptions import CplexError

from utility.Solution import Solution
from algorithms.Base import BiobjectiveSolver


class EpsilonGridWorker(mp.Process, BiobjectiveSolver):

    def __init__(self, z1, z2, inter_vars, task_q, done_q):
        mp.Process.__init__(self)
        BiobjectiveSolver.__init__(self, z1, z2, inter_vars)

        self.task_q = task_q
        self.done_q = done_q

    def run(self):

        while True:
            bound = self.task_q.get()
            if bound == "DONE":
                self.task_q.task_done()
                break

            try:
                z1_hat, r_hat = self._lexmin(0, 1, bound)
                self.done_q.put(z1_hat)
                self.task_q.task_done()
            except CplexError, exc:
                print exc
                self.task_q.task_done()
                continue

    def solve(self):
        raise NotImplementedError

    def approximate(self, gap):
        raise NotImplementedError