# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

'''
Asynchronous local execution.

Supports multicore architectures.
'''

from Bio.PopGen.Async import Async, FileRetriever

import thread

class Local(Async):
    '''Execution on Local machine.
    '''

    def __init__(self, num_cores = 1):
        '''Constructor.

           parameters:
           num_cores - Number of cores (for multiprocessor machines,
               multiply accordingly)
        '''
        Async.__init__(self)
        self.num_cores = num_cores
        self.cores_used = 0

    def _run_program(self, id, hook, parameters, input_files):
        '''Run program.

           For parameters, please check Async.run_program.

           Either runs a program if a core is available or
           schedules it.
        '''
        self.access_ds.acquire()
        self.waiting.append((id, hook, parameters, input_files))
        if self.cores_used < self.num_cores:
            self.cores_used += 1
            thread.start_new_thread(self.start_work, ())
        self.access_ds.release()

    def start_work(self):
        '''Starts work.

           Thread initial point.
           While there are tasks to be done, runs them.
           The thread dies as soon as there is nothing waiting to be
           executed.
        '''
        self.access_ds.acquire()
        while (len(self.waiting) > 0):
            id, hook, parameters, input_files = self.waiting[0]
            del self.waiting[0]
            self.running[id] = True
            self.access_ds.release()
            ret_code, output_files = hook.run_job(parameters, input_files)
            self.access_ds.acquire()
            del self.running[id]
            self.done[id] = ret_code, output_files
        self.cores_used -= 1
        self.access_ds.release()
  
