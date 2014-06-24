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
import os
from utility.Hypervolume import HyperVolume

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

        for i in xrange(2):
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

    def _terminate_worker(self, ids):
        """
        terminates the running worker on the cluster
        :return: None
        """
        cancel_call = self._config.get("CLUSTER", "cancel")
        commands = self._config.get("CLUSTER", "cancel_command")

        for id in ids:
            subprocess.call(cancel_call+" "+commands+" "+id)

    def _gently_terminate(self):
        """
            terminates the worker gently by sending an terminate signal
            to the worker
        """
        for _ in xrange(len(self.worker)):
            self.task_q.put_nowait("DONE")
        for _ in xrange(2):
            self.task_utopian_q.put(["DONE", None, None, None, None])
        time.sleep(3)
        self._manager.shutdown()

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

        self._gently_terminate()

        return self.solutions

    def approximate(self, gap, nof_sol=None):
        return self.solve(nof_sol=nof_sol)



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