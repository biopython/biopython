# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""
This module allows for the split (Async) of Fdist runs.
"""
import os
import thread
from time import sleep
from Bio.PopGen.Async import Local
from Bio.PopGen.FDist.Controller import FDistController

class FDistAsync(FDistController):
    def __init__(self, fdist_dir = '', ext = None):
        FDistController.__init__(self, fdist_dir, ext)

    def run_job(self, parameters, input_files):
        npops = parameters['npops']
        nsamples = parameters['nsamples']
        fst = parameters['fst']
        sample_size = parameters['sample_size']
        mut = parameters.get('mut', 0)
        num_sims = parameters.get('num_sims', 20000)
        data_dir = parameters.get('data_dir', '.')
        fst = self.run_fdist(npops, nsamples, fst, sample_size,
            mut, num_sims, data_dir)
        output_files = {}
        output_files['out.dat'] = open(data_dir + os.sep + 'out.dat', 'r')
        return fst, output_files

class SplitFDist:
    def __init__(self, report_fun = None,
        num_thr = 2, split_size = 1000, fdist_dir = '', ext = None):
        """

          report_fun - function to be called when a part of fdist is done
                       Important: Concurrent calls are possible
        """
        self.async = Local.Local(num_thr)
        self.async.hooks['fdist'] = FDistAsync(fdist_dir, ext)
        self.report_fun = report_fun
        self.split_size = split_size

    #There might be races when reporting...
    def monitor(self):
        while(True):
            sleep(1)
            self.async.access_ds.acquire()
            keys =  self.async.done.keys()[:]
            self.async.access_ds.release()
            for done in keys:
                self.async.access_ds.acquire()
                fst, files = self.async.done[done]
                del self.async.done[done]
                out_dat = files['out.dat']
                f = open(self.data_dir + os.sep + 'out.dat','a')
                f.writelines(out_dat.readlines())
                f.close()
                out_dat.close()
                self.async.access_ds.release()
                for file in os.listdir(self.parts[done]):
                    os.remove (self.parts[done] + os.sep + file)
                os.rmdir(self.parts[done])
                #print fst, out_dat
                if self.report_fun:
                    self.report_fun(fst)
            self.async.access_ds.acquire()
            if len(self.async.waiting) == 0 and len(self.async.running) == 0 \
               and len(self.async.done) == 0:
                break
            self.async.access_ds.release()
            #print 'R', self.async.running
            #print 'W', self.async.waiting
            #print 'R', self.async.running

    def acquire(self):
        self.async.access_ds.acquire()

    def release(self):
        self.async.access_ds.release()

    #You can only run a fdist case at a time
    def run_fdist(self, npops, nsamples, fst, sample_size,
        mut = 0, num_sims = 20000, data_dir='.'):
        num_parts = num_sims/self.split_size
        self.parts = {}
        self.data_dir = data_dir
        for directory in range(num_parts):
           full_path = data_dir + os.sep + str(directory)
           try:
               os.mkdir(full_path)
           except OSError:
               pass #Its ok, if it is already there
           id = self.async.run_program('fdist', {
               'npops'       : npops,
               'nsamples'    : nsamples,
               'fst'         : fst,
               'sample_size' : sample_size,
               'mut'         : mut,
               'num_sims'    : self.split_size,
               'data_dir'    : full_path
           }, {})
           self.parts[id] = full_path
        thread.start_new_thread(self.monitor, ())
