import itertools
import cplex
import numpy
import multiprocessing as mp

from utility.Solution import Solution


class RectangleSplittingWorker(mp.Process):

    def __init__(self, z1, z2, biob_cons, inter_vars, task_q, done_q, has_cons=False):
        mp.Process.__init__(self)
        self._models = (z1, z2)
        self._changeable_constraints = biob_cons
        self._inter_variables = inter_vars
        self._variables = z1.variables.get_names()
        self.task_q = task_q
        self.done_q = done_q

        #modify model for solving (absolutily inefficient but dont know how else)
        if has_cons:
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
            if z1_idx:
                sol, warm = self._lexmin(z1_idx, z2_idx, boundary,  warmstart=warmst[0], effort_level=0)
                self.done_q.put((z1_idx, sol, [warmst[0], warm], rectangle, hashs))
            else:
                sol, warm = self._lexmin(z1_idx, z2_idx, boundary,  warmstart=warmst[1], effort_level=0)
                self.done_q.put((z1_idx, sol, [warm, warmst[1]], rectangle, hashs))
            self.task_q.task_done()

