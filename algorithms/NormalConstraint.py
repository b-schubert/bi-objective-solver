import numpy
import cplex

from cplex.exceptions import CplexError

from Base import BiobjectiveSolver


class NormalConstraint(BiobjectiveSolver):
    """
        implements the normal constraint method from
        Messac et al., 2003, "The Normalized Normal Constraint Method for Generating the Pareto Frontier",
        Structural and Multidisciplinary Optimization
    """

    def __init__(self, z1, z2, inter_vars):
        super(NormalConstraint, self).__init__(z1, z2, inter_vars)
