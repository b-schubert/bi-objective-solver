__author__ = 'schubert'
from rectangle_worker import RectangleSplittingWorker
from multiprocessing.managers import SyncManager


class EpsilonRectangleWorker(RectangleSplittingWorker):

    def __init__(self, z1_name, z2_name, biob_cons, inter_vars, port, authkey, ip, nof_cpu=6, has_constraints=False):
        RectangleSplittingWorker.__init__(self, z1_name, z2_name, biob_cons, inter_vars, port, authkey, ip, nof_cpu=6, has_constraints=False)
        self.manager = self.__make_client_manager(port, authkey, ip)
        self.done_q = self.manager.get_done_q()
        self.task_q = self.manager.get_task_utopian_q()

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