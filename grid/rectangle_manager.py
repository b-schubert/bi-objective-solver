"""
    This file contains grid implementations of the rectangle splitting
    algorithm. The manager should be started on the master node of the grid
    and the worker are submitted as jobs and communicate via TCP/IP with the manager

    @Author: Benjamin Schubert

"""
from __future__ import division
import cplex
import heapq
import numpy
import argparse
import multiprocessing as mp
from multiprocessing.managers import SyncManager
import ConfigParser
import subprocess
import time

from algorithms.Base import BiobjectiveSolver

from utility.Hypervolume import HyperVolume


class RectangleSplittingManager(object):

    def __init__(self, z1_name, z2_name, inter_vars, nof_worker, nof_cpu=6, port=60, authkey="deimmuno", config="./config.cfg"):
        self.solutions = []
        self.models = (z1_name, z2_name)
        self.biob_cons = ["z2_cons", "z1_cons"]
        self.inter_vars = inter_vars
        self.nof_worker = nof_worker
        self.worker = []
        self._hypervol = HyperVolume()
        self._manager = self.__make_manager_server(port, authkey)
        #concurrent stuff
        self.task_q = self._manager.get_task_q()
        self.done_q = self._manager.get_done_q()
        self._config = ConfigParser.ConfigParser()
        self._config.read(config)
        self.empty_rectangles = []

        #setup worker
        self.__start_worker(nof_worker, port, authkey, nof_cpu)

    def __make_manager_server(self, port, authkey):
        """
         Starts a manager server
        :return: Manager object
        """

        job_q = mp.JoinableQueue()
        result_q = mp.Queue()

        # This is based on the examples in the official docs of multiprocessing.
        # get_{job|result}_q return synchronized proxies for the actual Queue
        # objects.
        class JobQueueManager(mp.managers.SyncManager):
            pass

        JobQueueManager.register('get_task_q', callable=lambda: job_q)
        JobQueueManager.register('get_done_q', callable=lambda: result_q)

        manager = JobQueueManager(address=('', port), authkey=authkey)
        manager.start()
        print 'Server started at port %s' % port
        return manager

    def __start_worker(self, nof_worker, port, authkey, nof_cpu):
        """
        Submit worker jobs to the queue.
        :param nof_cpu:
        :param port:
        :param authkey:
        :return: None
        """
        submit_call = self._config.get("CLUSTER", "submit")
        command = self._config.get("CLUSTER", "submit_command")

        template = open(self._config.get("GENERAL", "template"),"r").readlines()
        job_input_folder = self._config.get("GENERAL", "input")
        job_out_folder = self._config.get("GENERAL", "output")
        job_error_folder = job_out_folder = self._config.get("GENERAL", "error")
        rectangel_wroker = self._config.get("WORKER", "rectangle")

        for i in xrange(nof_worker):
            name = ".".join(self.models[0].split(".")[:-1])+"_worker_%i"%i
            input_log = job_input_folder+name+".sh"
            output_log = job_out_folder+name+".o"
            error_log = job_error_folder+name+".e"
            command = command%(nof_cpu, output_log, error_log, input_log)
            with open(input_log, "w") as f:
                f.write("".join(template))
                f.write("\n\npython %s -i %s %s -p %i -a %s -t %i -c %s -v %s"%(rectangel_wroker, self.models[0],
                                                                                 self.models[1], port, authkey,
                                                                                 nof_cpu, " ".join(self.biob_cons),
                                                                                 " ".join(self.inter_vars)))
            id = subprocess.check_output("%s %s "%(submit_call, command), shell=True)
            print id
            self.worker.append(id)

    def _terminate_worker(self):
        """
        terminates the running worker on the cluster
        :return: None
        """
        cancel_call = self._config.get("CLUSTER", "cancel")
        commands = self._config.get("CLUSTER", "cancel_command")

        for id in self.worker:
            subprocess.call(cancel_call+" "+commands+" "+id)

    def _gently_terminate(self):
        """
            terminates the worker gently by sending an terminate signal
            to the worker
        """
        for _ in xrange(len(self.worker)):
            self.task_q.put_nowait(["DONE", None, None, None, None])
        time.sleep(2)
        self._manager.shutdown()

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
            task_count = 2
            #init problems to solve
            self.task_q.put_nowait((0, 1, cplex.infinity, [None, None], ((None, None), (None, None))))
            self.task_q.put_nowait((1, 0, cplex.infinity, [None, None], ((None, None), (None, None))))

            self.task_q.join()


            b = [None, None]
            warm = [None, None]
            while task_count:
                pos, sol, warmstart, origin_rect = self.done_q.get()
                if pos is None:
                    continue
                task_count -= 1
                self.solutions.append(sol)
                b[pos] = sol.objs
                warm[pos] = warmstart[pos]

            rec_b = 0.5*(b[0][1]+b[1][1])
            task_count = 1
            self.task_q.put_nowait((0, 1, rec_b, warm, b))

        while task_count:

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
            task_count += 2
            #init problems to solve
            self.task_q.put_nowait((0, 1, cplex.infinity, [None, None], ((None, None), (None, None))))
            self.task_q.put_nowait((1, 0, cplex.infinity, [None, None], ((None, None), (None, None))))

            self.task_q.join()

            init_rect = [None, None]
            warm = [None, None]
            while task_count:
                pos, sol, warmstart, origin_rect = self.done_q.get()

                if pos is None:
                    continue
                task_count -= 1
                heapq.heappush(self.solutions, sol)
                init_rect[pos] = sol.objs
                warm[pos] = warmstart[pos]
            init_rect = tuple(init_rect)
            rec_b = 0.5*(init_rect[0][1]+init_rect[1][1])
            task_count += 1
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
                self._gently_terminate()
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

if __name__ == "__main__":
    import cPickle as pcl


    parser = argparse.ArgumentParser(description=' Rectangle Worker Grid implementation')
    parser.add_argument('--input','-i',
                      required=True,
                      nargs=2,
                      help="model files ")
    parser.add_argument('--output','-o',
                      required=True,
                      nargs=2,
                      help="Solution output as pickel")

    args = parser.parse_args()
    manager = RectangleSplittingManager(args.input[0],args.input[1],["x","y"],4)

    sols = manager.approximate(0.001)
    print sols
    pcl.dump(sols, open(args.output, "w"), -1)