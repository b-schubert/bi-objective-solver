import numpy
import itertools


class ParetoFilter(object):
    """
        This class implements a pareto filter
        which gets a list of solutions and filters for globally pareto optimal points

        See: Messac et al. 2003
    """
    @staticmethod
    def filter(solutions):
        """
            this function filters a list of solutions for globally pareto optimal points

            @param solutions: A list of solutions
            @type solutions: Solution
            @return: A list of global pareto optimal solutions
        """
        global_sols = []
        m = len(solutions)

        for i in xrange(m):
            zi = solutions[i]
            for j in xrange(m):
                if i == j:
                    continue

                zj = solutions[j]
                if not zi == zj and all(numpy.greater_equal(zi.objs, zj.objs)):
                    break

            global_sols.append(zi)
        return global_sols