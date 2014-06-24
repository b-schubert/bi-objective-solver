from __future__ import division
import cplex
import argparse
import numpy
import multiprocessing as mp
from multiprocessing.managers import SyncManager
import ConfigParser
import subprocess
import time
import os
from utility.Hypervolume import HyperVolume
from utility.ParetoFilter import ParetoFilter
from algorithms.Base import BiobjectiveSolver

class EpsilonGridManager(object):
    """
        This is prarallel implementation of
        the epsilon contraint method with precomputed boundaries
    """

    def __init__(self, z1_name, z2_name, inter_vars, nof_worker, nof_cpu=6, port=60, authkey="deimmuno", config="./config.cfg", constraints=None):
        """
            init function

            @param z1_name: The ILP of the first objective
            @param z2_name: The ILP of the second objective
            @param nof_cpu: THe number of Processes to start
        """

        self.solutions = []
        self.models = (z1_name, z2_name)
        self.biob_cons = ["z2_cons", "z1_cons"] if constraints is None else constraints
        self.inter_vars = inter_vars
        self.nof_worker = nof_worker
        self.worker = []
        self._hypervol = HyperVolume()
        self._manager = self.__make_manager_server(port, authkey)
        #concurrent stuff
        self.task_q = self._manager.get_task_q()
        self.done_q = self._manager.get_done_q()
        self.task_utopian_q = self._manager.get_task_utopian_q()

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
        task_utopian_q = mp.JoinableQueue()

        # This is based on the examples in the official docs of multiprocessing.
        # get_{job|result}_q return synchronized proxies for the actual Queue
        # objects.
        class JobQueueManager(mp.managers.SyncManager):
            pass

        JobQueueManager.register('get_task_q', callable=lambda: job_q)
        JobQueueManager.register('get_done_q', callable=lambda: result_q)
        JobQueueManager.register('get_task_utopian_q', callable=lambda: task_utopian_q)
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
        rectangel_wroker = self._config.get("WORKER", "rectangle_epsilon_grid")

        epsilon_worker = self._config.get("WORKER", "epsilon_grid")
        for i in xrange(nof_worker):
            name = ".".join(os.path.basename(self.models[0]).split(".")[:-1])+"_worker_%i"%i
            input_log = job_input_folder+name+".sh"
            output_log = job_out_folder+name+".o"
            error_log = job_error_folder+name+".e"
            com = command%(nof_cpu, error_log, output_log, input_log)
            with open(input_log, "w") as f:
                f.write("".join(template))
                f.write("\n\npython %s -i %s %s -p %i -a %s -t %i -c %s -v %s"%(epsilon_worker, self.models[0],
                                                                                 self.models[1], port, authkey,
                                                                                 nof_cpu, " ".join(self.biob_cons),
                                                                                 " ".join(self.inter_vars)))
            id = subprocess.check_output("%s %s "%(submit_call, com), shell=True)
            print id
            self.worker.append(id)


        for i in xrange(nof_worker-1):
            name = ".".join(os.path.basename(self.models[0]).split(".")[:-1])+"_worker_%i"%i
            input_log = job_input_folder+name+".sh"
            output_log = job_out_folder+name+".o"
            error_log = job_error_folder+name+".e"
            com = command%(nof_cpu, error_log, output_log, input_log)
            with open(input_log, "w") as f:
                f.write("".join(template))
                f.write("\n\npython %s -i %s %s -p %i -a %s -t %i -c %s -v %s"%(rectangel_wroker, self.models[0],
                                                                                 self.models[1], port, authkey,
                                                                                 nof_cpu, " ".join(self.biob_cons),
                                                                                 " ".join(self.inter_vars)))
            id = subprocess.check_output("%s %s "%(submit_call, com), shell=True)
            print id
            self.worker.append(id)


    def _terminate_worker(self, ids):
        """
        terminates the running worker on the cluster
        :return: None
        """
        cancel_call = self._config.get("CLUSTER", "cancel")
        commands = self._config.get("CLUSTER", "cancel_command")

        for id in ids:
            subprocess.call(cancel_call+" "+commands+" "+id)

    def _gently_terminate(self, grid=True):
        """
            terminates the worker gently by sending an terminate signal
            to the worker
        """
        if grid:
            for _ in xrange(len(self.worker)):
                self.task_q.put_nowait("DONE")
        else:
            for _ in xrange(len(self.worker)-1):
                self.task_utopian_q.put(["DONE", None, None, None, None])
        time.sleep(3)
        self._manager.shutdown()

    def solve_grid(self, nof_sol=None):


        solutions = []
        nof_sol = len(self.worker) if nof_sol is None else nof_sol
        #init problems to solve
        self.task_utopian_q.put((0, 1, cplex.infinity, [None, None], ((None, None), (None, None))))
        self.task_utopian_q.put((1, 0, cplex.infinity, [None, None], ((None, None), (None, None))))

        self.task_utopian_q.join()

        b = [None, None]
        warm = [None, None]
        while not self.done_q.empty():
            pos, sol, warmstart, origin_rect = self.done_q.get()

            solutions.append(sol)
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
            solutions.append(sols)

        self._gently_terminate(grid=True)
        return solutions

    def solve_rectangle(self, init_recs=None):
        """
            solve rectangle spliting after intial deterimnation fo rectangles
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
        self._gently_terminate(grid=False)
        self._terminate_worker()

    def solve(self, nof_sol=None):
        sols = self.solve_grid(nof_sol=nof_sol)

        #filter non pareto points
        sols = sorted(ParetoFilter.filter(sols))
        #self.solutions.extend(sols)

        print "Filtered pareto points"
        for s in sols:
            print s

        self.solve_rectangle(init_recs=sols)

        return self.solutions

    def approximate(self, gap, nof_sol=None):
        sols = self.solve_grid(nof_sol=nof_sol)

        #filter non pareto points
        sols = sorted(ParetoFilter.filter(sols))
        #self.solutions.extend(sols)

        print "Filtered pareto points"
        for s in sols:
            print s

        curr_gap = HyperVolume.calc_hypervol_gap(self.solutions, [sols[0].objs, sols[-1].objs], [])
        if numpy.allclose(gap, curr_gap, rtol=1e-03, atol=1e-04) or curr_gap < gap:
            return sols

        self.solve_rectangle(init_recs=sols)

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
    manager = EpsilonGridManager(args.input[0],args.input[1],["x","y"],4)

    sols = manager.approximate(0.001)
    print sols
    pcl.dump(sols, open(args.output, "w"), -1)