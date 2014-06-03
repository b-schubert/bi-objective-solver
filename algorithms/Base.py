import itertools
import cplex
import abc
import numpy

from utility.Solution import Solution
from utility.Hypervolume import HyperVolume


class BiobjectiveSolver(object):
    """
        Generic biobjectiv solver
    """

    __metaclass__ = abc.ABCMeta

    EPS = 1e-3
    #EPS = 1

    def __init__(self, z1, z2, interesting_vars, constraints=None):
        """
            initializes the solver model and
             modifies the given cplex model

             @param z1: cplex model with first objective
             @param z2: equivalent cplex model with second objective
             @param biob_constraints: a list of constraints which have to be changed in each iteration
                    should be only two!
             @param interesting_vars: list of variable names of the model which a user is interested in
        """
        z1 = cplex.Cplex(z1)
        z2 = cplex.Cplex(z2)
        z1.set_results_stream(None)
        z2.set_results_stream(None)

        self._models = (z1, z2)
        self._changeable_constraints = ["z2_cons", "z1_cons"] if constraints is None else constraints
        self._inter_variables = filter(lambda x:  x[0] in interesting_vars, z1.variables.get_names())
        self._solutions = []
        self._variables = z1.variables.get_names()
        self._hypervol = HyperVolume()

        #transform problems by adding an additional box constraint to each of the two
        #concerning the other objective function

        #dont know why i first have to go through the complete list and filter it
        if constraints is None:
            z1_obj_val = { v:z1.objective.get_linear(v) for v in self._variables if not numpy.allclose(z1.objective.get_linear(v), 0.0)}
            z2_obj_val = { v:z2.objective.get_linear(v) for v in self._variables if not numpy.allclose(z2.objective.get_linear(v),0.0)}
            z1.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=z2_obj_val.keys(), val=z2_obj_val.values())], senses=["L"], rhs=[0.0], range_values=[0], names=["z2_cons"])
            z2.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=z1_obj_val.keys(), val=z1_obj_val.values())], senses=["L"], rhs=[0.0], names=["z1_cons"])

    def _lexmin(self, z1_idx, z2_idx, boundary,  warmstart=None, effort_level=0):
        """
            lexicographic optimization based on input values

            @param z1_idx, z2_idx: integers defining which of the two objectives is optimized first
            @param z1_bound, boundary: floating point of new boundary of the
            @param warmstart (optional): is an optional feature to use a warmstart for z1
            @return Solution, values_of_all_variables_for_warm_start


            @QUESTION:Dont know if one should first delete old MIP starts or just keep adding them
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

    @abc.abstractmethod
    def solve(self):
        """
            solver functionality
        """
        NotImplementedError("solve function has to be implemented")

    @abc.abstractmethod
    def approximate(self, gap):
        """
            implements an approximation technique
        """
        NotImplementedError("solve function has to be implemented")

    def get_solutions(self):
        return self._solutions
