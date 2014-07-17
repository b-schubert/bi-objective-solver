from __future__ import division
import sys
import argparse
import ConfigParser
from multiprocessing.managers import SyncManager
from cplex.exceptions import CplexError
from Base import BiobjectiveSolver


class EpsilonGridWorker(BiobjectiveSolver):

    def __init__(self, z1, z2, biob_cons, inter_vars, port, authkey, ip, nof_cpu=6, has_constraints=False):
        BiobjectiveSolver.__init__(self, z1, z2, inter_vars, constraints=biob_cons if has_constraints else None,
                                   nof_cpu=nof_cpu)
        self.manager = self.__make_client_manager(port, authkey, ip)
        self.task_q = self.manager.get_task_q()
        self.done_q = self.manager.get_eps_done_q()

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
        ServerQueueManager.register('get_eps_done_q')

        manager = ServerQueueManager(address=(ip, port), authkey=authkey)
        manager.connect()

        print 'Client connected to %s:%s' % (ip, port)
        return manager

    def run(self):

        while True:
            print "waiting for task"
            bound = self.task_q.get()
            if bound == "DONE":
                self.task_q.task_done()
                break

            try:
                z1_hat, r_hat = self._lexmin(0, 1, bound)
                print z1_hat
                self.done_q.put(z1_hat)
                self.task_q.task_done()
            except CplexError, exc:
                print exc
                self.task_q.task_done()
                continue

    def solve(self):
        raise NotImplementedError

    def approximate(self, gap):
        raise NotImplementedError


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
                      type=int,
                      required=True,
                      help="nof of core"
                      )
    parser.add_argument('--constraints','-c',
                      required=True,
                      nargs=2,
                      help="Constraints"
                      )
    parser.add_argument('--variables','-v',
                      required=True,
                      nargs="+",
                      help="interesting variables"
                      )
    parser.add_argument('--hasconst','-hc',
                      action="store_true",
                      help="If constraints are already in model included"
                      )
    args = parser.parse_args()
    config = ConfigParser.ConfigParser()
    config.read("/home-link/zxmqy30/bi-objective-solver/grid/config.cfg")

    #z1_name, z2_name, biob_cons, inter_vars, port, authkey, ip, nof_cpu=6
    worker = EpsilonGridWorker(args.input[0], args.input[1],
                               args.constraints,args.variables, args.port, args.authkey,
                               config.get("GENERAL","master_ip"), args.threads, args.hasconst)

    sys.exit()
