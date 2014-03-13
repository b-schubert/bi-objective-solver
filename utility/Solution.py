import numpy


class Solution(object):
    """
        solution class for a mutliobjective problem
    """

    def __init__(self, objs, vars):
        """
            vars saves the variables which are of interest to the user and are not 0
        """
        self.objs = tuple(objs)
        self.vars = vars


    def __eq__(self, solution):
        """
            compare function with other solution
        """
        return numpy.allclose(self.objs, solution)

    def __lt__(self, solution):
        """
           compares a Solution based on its first obejctive value
        """
        return self.objs[0] <= solution.objs[0]

    def __str__(self):
        return str(self.objs)

    def __repr__(self):
        return str(self.objs)
