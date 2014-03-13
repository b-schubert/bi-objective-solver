from __future__ import division
import cplex
import numpy
import itertools

from copy import copy
from cplex.exceptions import CplexError


from Base import BiobjectiveSolver
from utility.Solution import Solution


class NormalConstraint(BiobjectiveSolver):
    """
        Implements the Normal Constraint Method to solve uniformly distributed pareto points
        for a bi-objective problem.

        See:
        (1) Messac et al., 2003, 'The normalized normal constraint method for generating the Pareto frontier',
            Struct Multidisc Optim 25, 86-98
    """

    def __init__(self, z1, z2, inter_vars):
        super(NormalConstraint, self).__init__(z1, z2, inter_vars)
        #anchor points


    def solve(self, nof_sol):
        """
            implements the algorithm.

            NOTICE: that the generated points do not have to be global pareto points
            if the problem is not convex!

            @param nof_sol: The number of solutions to generate
            @type nof_sol: int

        """

        def normalize(z):
            return numpy.array([(z[0] - z_t.objs[0])/l1, (z[1] - z_b.objs[1])/l2])

        def transform(z):
            return z[0]*l1+z_t.objs[0], z[1]*l2+z_b.objs[1]

        ##### setup ####
        #first generate anchor points
        z_t, r_t = self._lexmin(0, 1, cplex.infinity)
        z_b, r_b = self._lexmin(1, 0, cplex.infinity)



        self._solutions.append(z_t)
        self._solutions.append(z_b)

        #calculate l1, l2
        l1 = z_b.objs[0] - z_t.objs[0]
        l2 = z_t.objs[1] - z_b.objs[1]

        #normalize utopian points
        z_t_bar = normalize(z_t)
        z_b_bar = normalize(z_b)

        #generate utopian line vector
        N1 = z_b_bar - z_t_bar

        #generated normalized incremental step
        delta = 1/(nof_sol - 1)

        #generate utopian line points
        alphas = [delta*i for i in xrange(1, nof_sol - 1)]
        Xp = numpy.matrix([alphas[i]*z_t_bar - (1 - alphas[i])*z_b_bar for i in xrange(nof_sol - 2)])

        #prepare problem
        z1_hat_val = {v:self._models[0].objective.get_linear(v)/l1 for v in self._variables
                      if self._models[0].objective.get_linear(v) != 0}
        z2_hat_val = {v:self._models[1].objective.get_linear(v)/l1 for v in self._variables
                      if self._models[0].objective.get_linear(v) != 0}

        #dont know if one needs deep copy here
        z2 = copy(self._models[1])
        z2.objective.set_linear(z2_hat_val.iteritems())
        z2.linear_constraints.delete(self._changeable_constraints[1])
        z2.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=z1_hat_val.keys(),
                                                             val=numpy.multiply(z1_hat_val.values(), N1[0])),
                                            cplex.SparsePair(ind=z2_hat_val.keys(),
                                                             val=numpy.multiply(z2_hat_val.values(), N1[1]))],
                                  senses=["L", "L"],
                                  rhs=[cplex.infinity, cplex.infinity],
                                  range_values=[0.0, 0.0],
                                  names=["normal_cons_z1", "normal_cons_z2"])

        #solve dat shit!
        for i in xrange(nof_sol):
            z2.linear_constraints.set_rhs("normal_cons_z1", Xp[i, 0])
            z2.linear_constraints.set_rhs("normal_cons_z1", Xp[i, 1])
            try:
                r = z2.solve()

                obj = [numpy.inner(z1_hat_val.values(), r.solution.get_values(z1_hat_val.keys()))*l1,
                   r.solution.get_objective_value()*l2]

                inter_vars = {}
                if len(self._inter_variables) > 0:
                    inter_vars = {k:v for v, k in itertools.izip(r.solution.get_values(self._inter_variables),
                                                             self._inter_variables) if v > 0.0}
                self._solutions.append(Solution(obj, inter_vars))
            except CplexError, exc:
                print exc
                continue

        return self._solutions

    def approximate(self, gap, nof_solution):
        return self.solve(nof_solution)