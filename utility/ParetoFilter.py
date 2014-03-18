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
            is_global = True
            zi = solutions[i]
            for j in xrange(m):
                if i == j:
                    continue

                zj = solutions[j]
                print "Zi ", zi
                print "Zj ", zj
                print not zi == zj, all(numpy.greater_equal(zi.objs, zj.objs)), not zi == zj and all(numpy.greater_equal(zi.objs, zj.objs))
                if not zi == zj and all(numpy.greater_equal(zi.objs, zj.objs)):
                    is_global = False
                    break

            if is_global:
                global_sols.append(zi)
        return global_sols