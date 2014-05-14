"""
    This file contains concurrecnt implementations of the rectangle splitting
    algorithm

    @Author: Benjamin Schubert

"""

import cplex
import heapq
import numpy
import multiprocessing as mp

from algorithms.Base import BiobjectiveSolver
from rectangle_worker import RectangleSplittingWorker

from utility.Hypervolume import HyperVolume


class RectangleSplittingManager(object):

    def __init__(self, z1_name, z2_name, inter_vars, nof_worker):
        self.solutions = []
        self.models = (z1_name, z2_name)
        self.biob_cons = ["z2_cons", "z1_cons"]
        self.inter_vars = inter_vars
        self.nof_worker = nof_worker
        self.worker = []
        self._hypervol = HyperVolume()

        #concurrent stuff
        self.task_q = mp.JoinableQueue()
        self.done_q = mp.Queue()
        m = mp.Manager()

        self.empty_rectangles = m.list()
        #setup worker
        z_t = cplex.Cplex(z1_name)
        z_b = cplex.Cplex(z2_name)
        for _ in xrange(nof_worker):
            z1 = cplex.Cplex(z_t)
            z2 = cplex.Cplex(z_b)
            #z1.parameters.mip.strategy.search.set(1)
            #z2.parameters.mip.strategy.search.set(1)
            z1.parameters.threads.set(max(int(mp.cpu_count()/nof_worker), 1))
            z2.parameters.threads.set(max(int(mp.cpu_count()/nof_worker), 1))
            z1.set_results_stream(None)
            z2.set_results_stream(None)

            p = RectangleSplittingWorker(z1, z2, self.biob_cons, inter_vars, self.task_q, self.done_q)
            p.deamon = True
            self.worker.append(p)
            p.start()

    def _send_terminate(self):
        """
            terminates the worker processes prematurily
            only used in approximate to break the compuation
            after reaching a certain predefined gap
        """
        self.task_q.close()
        self.task_q.join_thread()
        self.done_q.close()
        self.done_q.join_thread()
        for p in self.worker:
            p.terminate()

    def _gently_terminate(self):
        """
            terminates the worker gently by sending an terminate signal
            to the worker
        """
        for _ in xrange(len(self.worker)):
            self.task_q.put_nowait(["DONE", None, None, None, None])
        self.task_q.close()
        self.task_q.join()
        for p in self.worker:
            print "shut down worker ", p.pid
            p.join()

    def solve(self, init_recs=None):
        """
            coordinates the solving step
        """

        task_count = 0
        if init_recs:
            b = [init_recs[0].objs, init_recs[-1].objs]
            self.solutions.extend(init_recs)
            for i in xrange(len(init_recs)-1):
                task_count += 1
                zi = init_recs[i]
                zj = init_recs[i+1]
                warm = [zi.warm_start, zj.warm_start]
                rec = [zi.objs, zj.objs]
                rec_b = 0.5*(rec[0][1]+rec[1][1])
                self.task_q.put_nowait((0, 1, rec_b, warm, rec))

        else:
            task_count = 1
            #init problems to solve
            self.task_q.put_nowait((0, 1, cplex.infinity, [None, None], ((None, None), (None, None))))
            self.task_q.put_nowait((1, 0, cplex.infinity, [None, None], ((None, None), (None, None))))

            self.task_q.join()


            b = [None, None]
            warm = [None, None]
            while not self.done_q.empty():
                pos, sol, warmstart, origin_rect = self.done_q.get()
                self.solutions.append(sol)
                b[pos] = sol.objs
                warm[pos] = warmstart[pos]

            rec_b = 0.5*(b[0][1]+b[1][1])

            self.task_q.put_nowait((0, 1, rec_b, warm, b))

        while task_count > 0:

            pos, sol, warm, origin_rect = self.done_q.get()
            print "Current Rectangle ", origin_rect
            print "Solution ", sol
            print "Solutions ", self.solutions

            #lexmin2
            if pos:
                print "lexmin2", pos
                if not numpy.allclose(sol.objs, origin_rect[0], rtol=1e-01, atol=1e-04):
                    rec = (origin_rect[0], sol.objs)
                    self.solutions.append(sol)
                    rec_b = 0.5*(rec[0][1]+rec[1][1])
                    self.task_q.put_nowait((0, 1, rec_b, warm, rec))
                    task_count += 1
            #lexmin1
            else:
                print "lexmin1 ",pos
                rec_t = sol.objs[0]-BiobjectiveSolver.EPS
                self.task_q.put_nowait((1, 0, rec_t, warm, origin_rect))
                task_count += 1
                if not numpy.allclose(sol.objs, origin_rect[1], rtol=1e-01, atol=1e-04):
                    rec = (sol.objs, origin_rect[1])
                    rec_b = 0.5*(rec[0][1]+rec[1][1])
                    self.solutions.append(sol)
                    self.task_q.put_nowait((0, 1, rec_b, warm, rec))
                    task_count += 1

            task_count -= 1
            print "Tasks still running: ", task_count

        #all work done send exit signal
        self._gently_terminate()

        return self.solutions

    def approximate(self, gap, init_recs=None):
        task_count = 0
        proof_not_empty_rec = {}
        empty_rect = set()

        if init_recs:
            init_rect = [init_recs[0].objs, init_recs[-1].objs]
            self.solutions.extend(init_recs)
            for i in xrange(len(init_recs)-1):
                task_count += 1
                zi = init_recs[i]
                zj = init_recs[i+1]
                warm = [zi.warm_start, zj.warm_start]
                rec = [zi.objs, zj.objs]
                rec_b = 0.5*(rec[0][1]+rec[1][1])
                self.task_q.put_nowait((0, 1, rec_b, warm, rec))
        else:
            task_count += 1
            #init problems to solve
            self.task_q.put_nowait((0, 1, cplex.infinity, [None, None], ((None, None), (None, None))))
            self.task_q.put_nowait((1, 0, cplex.infinity, [None, None], ((None, None), (None, None))))

            self.task_q.join()

            init_rect = [None, None]
            warm = [None, None]
            while not self.done_q.empty():
                pos, sol, warmstart, origin_rect = self.done_q.get_nowait()
                heapq.heappush(self.solutions, sol)
                init_rect[pos] = sol.objs
                warm[pos] = warmstart[pos]
            init_rect = tuple(init_rect)
            rec_b = 0.5*(init_rect[0][1]+init_rect[1][1])

            self.task_q.put_nowait((0, 1, rec_b, warm, init_rect))

        while task_count > 0:

            pos, sol, warm, origin_rect = self.done_q.get()
            print
            print
            print "Current Rectangle ", origin_rect
            print "Solution: ", sol
            #print "Solutions ", self.solutions

            #lexmin2
            if pos:
                #print "lexmin2", pos
                if not numpy.allclose(sol.objs, origin_rect[0], rtol=1e-01, atol=1e-04):
                    rec = (origin_rect[0], sol.objs)
                    heapq.heappush(self.solutions, sol)
                    rec_b = 0.5*(rec[0][1]+rec[1][1])
                    self.task_q.put_nowait((0, 1, rec_b, warm, rec))
                    task_count += 1

                    if origin_rect in proof_not_empty_rec:
                        empty_rect.add((sol.objs, proof_not_empty_rec[origin_rect]))
                        del proof_not_empty_rec[origin_rect]
                    else:
                        #empty_rect.add((sol.objs, origin_rect[1]))
                        proof_not_empty_rec[origin_rect] = sol.objs

                else:
                    if origin_rect in proof_not_empty_rec:
                        empty_rect.add((origin_rect[0], proof_not_empty_rec[origin_rect]))
                        del proof_not_empty_rec[origin_rect]
                    else:
                        proof_not_empty_rec[origin_rect] = origin_rect[0]

            #lexmin1
            else:
                #print "lexmin1 ",pos
                rec_t = sol.objs[0]-BiobjectiveSolver.EPS
                self.task_q.put_nowait((1, 0, rec_t, warm, origin_rect))
                task_count += 1

                if not numpy.allclose(sol.objs, origin_rect[1], rtol=1e-01, atol=1e-04):
                    rec = (sol.objs, origin_rect[1])
                    rec_b = 0.5*(rec[0][1]+rec[1][1])
                    heapq.heappush(self.solutions, sol)
                    self.task_q.put_nowait((0, 1, rec_b, warm, rec))
                    task_count += 1
                    if origin_rect in proof_not_empty_rec:
                        empty_rect.add((proof_not_empty_rec[origin_rect], sol.objs))
                        del proof_not_empty_rec[origin_rect]

                    else:
                        proof_not_empty_rec[origin_rect] = sol.objs
                else:
                    if origin_rect in proof_not_empty_rec:
                        empty_rect.add((proof_not_empty_rec[origin_rect], origin_rect[1]))
                        del proof_not_empty_rec[origin_rect]

                    else:
                        proof_not_empty_rec[origin_rect] = origin_rect[1]

            task_count -= 1
            print "Tasks still running: ", task_count

            #if desired approximation quality is reached terminate and return results
            #print "Empty Rectangle ", empty_rect
            #print "Proofs ", proof_not_empty_rec
            cur_gap = self._hypervol.calc_hypervol_gap(self.solutions, init_rect, empty_rect)
            print "Current Hypervol gap is: ", cur_gap
            if numpy.allclose(gap, cur_gap, rtol=1e-01, atol=1e-04) or cur_gap < gap:
                self._send_terminate()
                return self.solutions
        print
        print
        print "Last hypervol gap: ", self._hypervol.calc_hypervol_gap(self.solutions, init_rect, empty_rect)
        #print "Empty rectangles: ", empty_rect
        #actually all worke is done terminate normally
        for _ in xrange(len(self.worker)):
            self.task_q.put_nowait(["DONE", None, None, None, None])
        self.task_q.close()
        self.task_q.join()
        for p in self.worker:
            print "shut down worker ", p.pid
            p.join()

        return self.solutions


class RectangleSplittingManagerDirectSplit(RectangleSplittingManager):
    '''
        a version of the splitter in which it is not waiting for lexmin_z1 to finish before lexmin_z2 starts
    '''

    def solve(self):
        task_count = 2

        #init problems to solve
        self.task_q.put_nowait((0, 1, cplex.infinity, None, ((None, None), (None, None))))
        self.task_q.put_nowait((1, 0, cplex.infinity, None, ((None, None), (None, None))))

        self.task_q.join()


        b = [None, None]
        warm = [None, None]
        while not self.done_q.empty():
            pos, sol, warmstart, origin_rect = self.done_q.get()
            self.solutions.append(sol)
            b[pos] = sol.objs
            warm[pos] = warmstart

        rec_b = 0.5*(b[0][1]+b[1][1])
        rec_t = b[1][0]-BiobjectiveSolver.EPS
        self.task_q.put_nowait((0, 1, rec_b, warm[1], b))
        self.task_q.put_nowait((1,0,rec_t, warm[0], b))

        while task_count > 0:
            pos, sol, warmstart, origin_rect = self.done_q.get()
            print "Current Rectangle ", origin_rect
            print "Solution: ", sol
            #print "Solutions ", self.solutions

            #lexmin2
            if pos:
                print "lexmin2", pos
                if not numpy.allclose(sol.objs, origin_rect[0]):
                    rec = (origin_rect[0], sol.objs)
                    self.solutions.append(sol)
                    rec_b = 0.5*(rec[0][1]+rec[1][1])
                    rec_t = rec[1][0]-BiobjectiveSolver.EPS
                    self.task_q.put_nowait((0, 1, rec_b, warmstart, rec))
                    self.task_q.put_nowait((1,0,rec_t, warmstart, rec))
                    task_count +=2
            else:
                if not numpy.allclose(sol.objs, origin_rect[1]):
                    rec = (sol.objs, origin_rect[1])
                    rec_b = 0.5*(rec[0][1]+rec[1][1])
                    rec_t = rec[1][0]-BiobjectiveSolver.EPS
                    self.solutions.append(sol)
                    self.task_q.put_nowait((0, 1, rec_b, warmstart, rec))
                    self.task_q.put_nowait((1,0,rec_t, warmstart, rec))
                    task_count += 2

            task_count -= 1
            print "Tasks still running: ", task_count

        #all work done send exit signal
        for _ in xrange(len(self.worker)):
            self.task_q.put_nowait(["DONE", None, None, None, None])
        self.task_q.close()
        self.task_q.join()
        for p in self.worker:
            print "shut down worker ", p.pid
            p.join()

        return self.solutions