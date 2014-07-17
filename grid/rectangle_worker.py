from __future__ import division
import itertools
import sys
import cplex
import numpy
import argparse
import ConfigParser
from multiprocessing.managers import SyncManager
from Solution import Solution

from cplex.exceptions import CplexError

class RectangleSplittingWorker(object):

    def __init__(self, z1_name, z2_name, biob_cons, inter_vars, port, authkey, ip, nof_cpu=6):
        z1 = cplex.Cplex(z1_name)
        z2 = cplex.Cplex(z2_name)
        print nof_cpu
        z1.parameters.threads.set(int(nof_cpu))
        z2.parameters.threads.set(int(nof_cpu))
        #z1.set_results_stream(None)
        #z2.set_results_stream(None)

        self.manager = self.__make_client_manager(port, authkey, ip)
        self._models = (z1, z2)
        self._changeable_constraints = biob_cons
        self._inter_variables = filter(lambda x:  x[0] in inter_vars, z1.variables.get_names())
        self._variables = z1.variables.get_names()
        self.task_q = self.manager.get_task_q()
        self.done_q = self.manager.get_done_q()

        #modify model for solving (absolutily inefficient but dont know how else)
        z1_obj_val = {}
        z2_obj_val = {}
        for v in self._variables:
            val_z1 = z1.objective.get_linear(v)
            val_z2 = z2.objective.get_linear(v)
            if not numpy.allclose(val_z1, 0.0):
                z1_obj_val[v] = val_z1
            if not numpy.allclose(val_z2, 0.0):
                z2_obj_val[v] = val_z2
        z1.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=z2_obj_val.keys(), val=z2_obj_val.values())], senses=["L"], rhs=[0.0], range_values=[0], names=[biob_cons[0]])
        z2.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=z1_obj_val.keys(), val=z1_obj_val.values())], senses=["L"], rhs=[0.0], names=[biob_cons[1]])

        #run the worker
        self.run()

    def __make_client_manager(self, port, authkey, ip):
        """
        generates manager and connects to manager master manager
        :param port:
        :param authkey:
        :return: manager
        """
        class ServerQueueManager(SyncManager):
            pass

        ServerQueueManager.register('get_task_q')
        ServerQueueManager.register('get_done_q')

        manager = ServerQueueManager(address=(ip, port), authkey=authkey)
        manager.connect()

        print 'Client connected to %s:%s' % (ip, port)
        return manager

    def _lexmin(self, z1_idx, z2_idx, boundary,  warmstart=None, effort_level=0):
        """
            lexicographic optimization based on input values

            :param z1_idx, z2_idx: integers defining which of the two objectives is optimized first
            :param z1_bound, boundary: floating point of new boundary of the
            :param warmstart (optional): is an optional feature to use a warmstart for z1
            :return Solution, values_of_all_variables_for_warm_start


            :QUESTION:Dont know if one should first delete old MIP starts or just keep adding them
        """

        z1 = self._models[z1_idx]
        if warmstart:
            z1.MIP_starts.delete() #this is questionable
            z1.MIP_starts.add([self._variables, warmstart], effort_level)
        z1.linear_constraints.set_rhs(self._changeable_constraints[z1_idx], boundary)
        z1.solve()
        z1_hat = z1.solution.get_objective_value()
        z1_hat_values = z1.solution.get_values()

        #lexicographical solution second model
        z2 = self._models[z2_idx]

        z2.linear_constraints.set_rhs(self._changeable_constraints[z2_idx], z1_hat)
        #set warm start from first objective
        z2.MIP_starts.add([self._variables, z1_hat_values], z2.MIP_starts.effort_level.auto)
        z2.solve()
        z2_hat = z2.solution.get_objective_value()
        z2_hat_values = z2.solution.get_values()

        #generate solution object
        objs = [0, 0]
        objs[z1_idx] = z1_hat
        objs[z2_idx] = z2_hat

        inter_vars = {}
        if len(self._inter_variables) > 0:
            inter_vars = {k:v for v, k in itertools.izip(z2.solution.get_values(self._inter_variables),
                                                         self._inter_variables)
                          if numpy.greater(v, 0.0) and not numpy.allclose(v, 0.0, rtol=1e-01, atol=1e-04)}
        s = Solution(objs, inter_vars)

        return s, z2_hat_values

    def run(self):
        while True:
            z1_idx, z2_idx, boundary, warmst, rectangle, hashs = self.task_q.get()
            if z1_idx == "DONE":
                self.task_q.task_done()
                break
            try:
                if z1_idx:
                        sol, warm = self._lexmin(z1_idx, z2_idx, boundary,  warmstart=warmst[0], effort_level=0)
                        self.done_q.put((z1_idx, sol, [warmst[0], warm], rectangle, hashs))
                else:
                        sol, warm = self._lexmin(z1_idx, z2_idx, boundary,  warmstart=warmst[1], effort_level=0)
                        self.done_q.put((z1_idx, sol, [warm, warmst[1]], rectangle,hashs))

                self.task_q.task_done()
            except CplexError, e:
                self.done_q.put(("error", None, None, None, None))
                self.task_q.task_done()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' Recnagle Worker Grid implementation')
    parser.add_argument('--input','-i',
                      required=True,
                      nargs=2,
                      help="model files ")
    parser.add_argument('--port','-p',
                      type=int,
                      required=True,
                      help="port to connect"
                      )
    parser.add_argument('--authkey','-a',
                      required=True,
                      help="authentication key"
                      )
    parser.add_argument('--threads','-t',
                      required=True,
                      type=int,
                      help="nof of core"
                      )
    parser.add_argument('--constraints','-c',
                      nargs=2,
                      required=True,
                      help="Constraints"
                      )
    parser.add_argument('--variables','-v',
                      nargs=2,
                      required=True,
                      help="interesting variables"
                      )
    args = parser.parse_args()
    config = ConfigParser.ConfigParser()
    config.read("/share/usr/schubert/projects/multiobjective_optimization_coopr/src/bi-objective-solver/bi-objective-solver/grid/config.cfg")

    #z1_name, z2_name, biob_cons, inter_vars, port, authkey, ip, nof_cpu=6
    worker = RectangleSplittingWorker(args.input[0], args.input[1],
                                      args.constraints,args.variables, args.port, args.authkey,
                                      config.get("GENERAL","master_ip"), args.threads)

    sys.exit()


