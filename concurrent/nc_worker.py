__author__ = 'schubert'
import numpy
import cplex
import itertools
import multiprocessing as mp

from cplex.exceptions import CplexError

from utility.Solution import Solution


class NormalConstraintWorker(mp.Process):

    def __init__(self, z1, z2, inter_vars, task_q, done_q):
        mp.Process.__init__(self)
        self._models = (z1, z2)
        self._variables = z1.variables.get_names()
        self._inter_variables = filter(lambda x:  x[0] in inter_vars, self._variables)
        self.task_q = task_q
        self.done_q = done_q

        #modify model for solving (absolutily inefficient but dont know how else)
        self._z1_obj_val = {v:z1.objective.get_linear(v) for v in self._variables if z1.objective.get_linear(v) != 0.0}
        self._z2_obj_val = {v:z2.objective.get_linear(v) for v in self._variables if z2.objective.get_linear(v) != 0.0}

    def run(self):
        modified = False
        while True:
            Xp, N1, warm = self.task_q.get()
            if Xp == "DONE":
                self.task_q.task_done()
                break

            z2 = self._models[1]
            z2.MIP_starts.delete()
            if warm:
                z2.MIP_starts.delete() #this is questionable
                z2.MIP_starts.add([self._variables, warm], 0)

            if not modified:
                z2.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=self._z1_obj_val.keys(),
                                                        val=numpy.multiply(self._z1_obj_val.values(), N1[0])),
                                                    cplex.SparsePair(ind=self._z2_obj_val.keys(),
                                                        val=numpy.multiply(self._z2_obj_val.values(), N1[1]))],
                                          senses=["L", "L"],
                                          rhs=[cplex.infinity, cplex.infinity],
                                          range_values=[0.0, 0.0],
                                          names=["normal_cons_z1", "normal_cons_z2"])
                modified = True

            z2.linear_constraints.set_rhs("normal_cons_z1", N1[0]*Xp[0, 0])
            z2.linear_constraints.set_rhs("normal_cons_z2", N1[1]*Xp[0, 1])

            try:
                z2.solve()
                #print "solved"
                obj = [numpy.inner(self._z1_obj_val.values(), z2.solution.get_values(self._z1_obj_val.keys())),
                       z2.solution.get_objective_value()]
                #print "Solution ", obj
                inter_vars = {}
                if len(self._inter_variables) > 0:
                    inter_vars = {k:v for v, k in itertools.izip(z2.solution.get_values(self._inter_variables),
                                                                 self._inter_variables) if v > 0.0}
                self.done_q.put(Solution(obj, inter_vars, warmstart=z2.solution.get_values()))
                self.task_q.task_done()
            except CplexError, exc:
                print exc
                self.task_q.task_done()
                continue