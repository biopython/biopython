# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Asynchronous execution of Fdist and spliting of loads.

FDistAsync Allows for the execution of FDist.

SplitFDist splits a single Fdist execution in several, taking advantage
of multi-core architectures.
"""

import os
import shutil
import threading
from time import sleep
from Bio.PopGen.Async import Local
from Bio.PopGen.FDist.Controller import FDistController

__docformat__ = "restructuredtext en"

class FDistAsync(FDistController):
    """Asynchronous FDist execution.
    """

    def __init__(self, fdist_dir="", ext=None):
        """Constructor.

        Parameters:
        
          - fdist_dir - Where fdist can be found, if = "", then it
              should be on the path.
          - ext - Extension of binary names (e.g. nothing on Unix,
                ".exe" on Windows
        """
        FDistController.__init__(self, fdist_dir, ext)

    def run_job(self, parameters, input_files):
        """Runs FDist asynchronously.

           Gets typical Fdist parameters from a dictionary and
           makes a "normal" call. This is run, normally, inside
           a separate thread.
        """
        npops = parameters['npops']
        nsamples = parameters['nsamples']
        fst = parameters['fst']
        sample_size = parameters['sample_size']
        mut = parameters.get('mut', 0)
        num_sims = parameters.get('num_sims', 20000)
        data_dir = parameters.get('data_dir', '.')
        is_dominant = parameters.get('is_dominant', False)
        theta = parameters.get('theta', 0.06)
        beta = parameters.get('beta', (0.25, 0.25))
        max_freq = parameters.get('max_freq', 0.99)
        fst = self.run_fdist(npops, nsamples, fst, sample_size,
                             mut, num_sims, data_dir,
                             is_dominant, theta, beta,
                             max_freq)
        output_files = {}
        output_files['out.dat'] = open(data_dir + os.sep + 'out.dat', 'r')
        return fst, output_files


class SplitFDist(object):
    """Splits a FDist run.

       The idea is to split a certain number of simulations in smaller
       numbers (e.g. 30.000 sims split in 30 packets of 1.000). This
       allows to run simulations in parallel, thus taking advantage
       of multi-core CPUs.

       Each SplitFDist object can only be used to run a single FDist
       simulation.
    """
    def __init__(self, report_fun=None,
                 num_thr=2, split_size=1000, fdist_dir='', ext=None):
        """Constructor.

           Parameters:
           
             - report_fun - Function that is called when a single packet is
               run, it should have a single parameter: Fst.
             - num_thr - Number of desired threads, typically the number
               of cores.
             - split_size - Size that a full simulation will be split in.
             - ext - Binary extension name (e.g. nothing on Unix, '.exe' on
               Windows).
        """
        self.async = Local.Local(num_thr)
        self.async.hooks['fdist'] = FDistAsync(fdist_dir, ext)
        self.report_fun = report_fun
        self.split_size = split_size

    # There might be races when reporting...
    def monitor(self):
        """Monitors and reports (using report_fun) execution.

           Every time a partial simulation ends, calls report_fun.
           IMPORTANT: monitor calls can be concurrent with other
           events, ie, a tasks might end while report_fun is being
           called. This means that report_fun should be consider that
           other events might be happening while it is running (it
           can call acquire/release if necessary).
        """
        while(True):
            sleep(1)
            self.async.access_ds.acquire()
            keys = list(self.async.done.keys())  # copy it
            self.async.access_ds.release()
            for done in keys:
                self.async.access_ds.acquire()
                fst, files = self.async.done[done]
                del self.async.done[done]
                out_dat = files['out.dat']
                f = open(self.data_dir + os.sep + 'out.dat', 'a')
                f.writelines(out_dat.readlines())
                f.close()
                out_dat.close()
                self.async.access_ds.release()
                for file in os.listdir(self.parts[done]):
                    os.remove(self.parts[done] + os.sep + file)
                os.rmdir(self.parts[done])
                if self.report_fun:
                    self.report_fun(fst)
            self.async.access_ds.acquire()
            if len(self.async.waiting) == 0 and len(self.async.running) == 0 \
               and len(self.async.done) == 0:
                break
            self.async.access_ds.release()

    def acquire(self):
        """Allows the external acquisition of the lock.
        """
        self.async.access_ds.acquire()

    def release(self):
        """Allows the external release of the lock.
        """
        self.async.access_ds.release()

    # You can only run a fdist case at a time
    def run_fdist(self, npops, nsamples, fst, sample_size,
                  mut=0, num_sims=20000, data_dir='.',
                  is_dominant=False, theta=0.06, beta=(0.25, 0.25),
                  max_freq=0.99):
        """Runs FDist.

           Parameters can be seen on FDistController.run_fdist.

           It will split a single execution in several parts and
           create separated data directories.
        """
        num_parts = num_sims // self.split_size
        self.parts = {}
        self.data_dir = data_dir
        for directory in range(num_parts):
            full_path = data_dir + os.sep + str(directory)
            try:
                os.mkdir(full_path)
            except OSError:
                pass  # Its ok, if it is already there
            if "ss_file" in os.listdir(data_dir):
                shutil.copy(data_dir + os.sep + "ss_file", full_path)
            id = self.async.run_program('fdist', {
                'npops': npops,
                'nsamples': nsamples,
                'fst': fst,
                'sample_size': sample_size,
                'mut': mut,
                'num_sims': self.split_size,
                'data_dir': full_path,
                'is_dominant': is_dominant,
                'theta': theta,
                'beta': beta,
                'max_freq': max_freq
            }, {})
            self.parts[id] = full_path
        threading.Thread(target=self.monitor).run()
