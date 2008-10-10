# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.



"""
This module allows to control fdist.

This will allow to call fdist and associated program (cplot, datacal, pv).

http://www.rubic.rdg.ac.uk/~mab/software.html
"""

import os
import tempfile
from sys import platform, maxint
from shutil import copyfile
from random import randint, random
from time import strftime, clock
#from logging import debug

class FDistController:
    def __init__(self, fdist_dir = '', ext = None):
        """Initializes the controller.
        
        fdist_dir is the directory where fdist2 is.
        ext is the extension of binaries (.exe on windows, 
          none on Unix)
        
        """
        self.tmp_idx = 0
        self.fdist_dir = fdist_dir
        self.os_name = os.name
        if self.os_name=='nt' or platform=='cygwin':
            py_ext = '.exe'
        else:
            py_ext = ''
        if ext == None:
            self.ext = py_ext
        else:
            self.ext = ext
        exec_counts = 0

    def _get_path(self, app):
        """Returns the path to an fdist application.

           Includes Path where fdist can be found plus executable extension.
        """
        if self.fdist_dir == '':
            return app + self.ext
        else:
            return os.sep.join([self.fdist_dir, app]) + self.ext

    def _get_temp_file(self):
        """Gets a temporary file name.

           Returns a temporary file name, if executing inside jython
           tries to replace unexisting tempfile.mkstemp().
        """
        self.tmp_idx += 1
        return strftime("%H%M%S") + str(int(clock()*100)) + str(randint(0,1000)) + str(self.tmp_idx)

    def run_datacal(self, data_dir='.'):
        """Executes datacal.
        
           data_dir - Where the data is found.
        """
        in_name = self._get_temp_file()
        out_name = self._get_temp_file()
        f = open(data_dir + os.sep + in_name, 'w')
        f.write('a\n')
        f.close()
        curr_dir = os.getcwd()
        os.system('cd ' + data_dir + ' && ' +
                self._get_path('datacal') + ' < ' + in_name + ' > ' + out_name)
        f = open(data_dir + os.sep + out_name)
        fst_line = f.readline().rstrip().split(' ')
        fst = float(fst_line[4])
        sample_line = f.readline().rstrip().split(' ')
        sample = int(sample_line[9])
        f.close()
        os.remove(data_dir + os.sep + in_name)
        os.remove(data_dir + os.sep + out_name)
        return fst, sample

    def _generate_intfile(self, data_dir):
        """Generates an INTFILE.

           Parameter:
           data_dir - data directory
        """
        inf = open(data_dir + os.sep + 'INTFILE', 'w')
        for i in range(98):
            inf.write(str(randint(-maxint+1,maxint-1)) + '\n') 
        inf.write('8\n')
        inf.close()
    
    def run_fdist(self, npops, nsamples, fst, sample_size,
        mut = 0, num_sims = 20000, data_dir='.'):
        """Executes fdist.
        
        Parameters:
        npops - Number of populations
        nsamples - Number of populations sampled
        fst - expected Fst
        sample_size - Sample size per population
        mut - 1=Stepwise, 0=Infinite allele
        num_sims - number of simulations
        data_dir - Where the data is found

        Returns:
        fst - Average Fst
        
        Important Note: This can take quite a while to run!
        """
        if fst >= 0.9:
            #Lets not joke
            fst = 0.899
        if fst <= 0.0:
            #0  will make fdist run forever
            fst = 0.001
        in_name = 'input.fd'
        out_name = 'output.fd'
        #print 'writing', data_dir + os.sep + in_name
        f = open(data_dir + os.sep + in_name, 'w')
        f.write('y\n\n')
        f.close()
        f = open(data_dir + os.sep + 'fdist_params2.dat', 'w')
        f.write(str(npops) + '\n')
        f.write(str(nsamples) + '\n')
        f.write(str(fst) + '\n')
        f.write(str(sample_size) + '\n')
        f.write(str(mut) + '\n')
        f.write(str(num_sims) + '\n')
        f.close()
        self._generate_intfile(data_dir)

        os.system('cd ' + data_dir + ' && ' +
            self._get_path('fdist2') + ' < ' + in_name + ' > ' + out_name)
        f = open(data_dir + os.sep + out_name)
        lines = f.readlines()
        f.close()
        for line in lines:
          if line.startswith('average Fst'):
            fst = float(line.rstrip().split(' ')[-1])
        os.remove(data_dir + os.sep + in_name)
        os.remove(data_dir + os.sep + out_name)
        return fst

    def run_fdist_force_fst(self, npops, nsamples, fst, sample_size,
        mut = 0, num_sims = 20000, data_dir='.', try_runs = 5000, limit=0.001):
        """Exectues fdist trying to force Fst.
        
        Parameters:
        try_runs - Number of simulations on the part trying to get
                   Fst correct
        limit - Interval limit
        Other parameters can be seen on run_fdist.
        """
        max_run_fst = 1
        min_run_fst = 0
        current_run_fst = fst
        old_fst = fst
        while True:
            #debug('testing fst ' +  str(current_run_fst))
            real_fst = self.run_fdist(npops, nsamples, current_run_fst, sample_size,
                mut, try_runs, data_dir)
            #debug('got real fst ' +  str(real_fst))
            if abs(real_fst - fst) < limit:
                #debug('We are OK')
                return self.run_fdist(npops, nsamples, current_run_fst, sample_size,
                    mut, num_sims, data_dir)
            old_fst = current_run_fst
            if real_fst > fst:
                max_run_fst = current_run_fst
                if current_run_fst < min_run_fst + limit:
                    #we can do no better
                    #debug('Lower limit is ' + str(min_run_fst))
                    return self.run_fdist(npops, nsamples, current_run_fst,
                        sample_size, mut, num_sims, data_dir)
                current_run_fst = (min_run_fst + current_run_fst)/2
            else:
                min_run_fst = current_run_fst
                if current_run_fst > max_run_fst - limit:
                    #we can do no better
                    #debug('Upper limit is ' + str(max_run_fst))
                    return self.run_fdist(npops, nsamples, current_run_fst,
                        sample_size, mut, num_sims, data_dir)
                current_run_fst = (max_run_fst + current_run_fst)/2

    def run_cplot(self, ci= 0.95, data_dir='.'):
        """Executes cplot.

        ci - Confidence interval.
        data_dir - Where the data is found.
        """
        in_name = self._get_temp_file()
        out_name = self._get_temp_file()
        f = open(data_dir + os.sep + in_name, 'w')
        f.write('out.dat out.cpl\n' + str(ci) + '\n')
        f.close()
        curr_dir = os.getcwd()
        self._generate_intfile(data_dir)
        os.system('cd ' + data_dir + ' && '  +
            self._get_path('cplot') + ' < ' + in_name + ' > ' + out_name)
        os.remove(data_dir + os.sep + in_name)
        os.remove(data_dir + os.sep + out_name)
        f = open(data_dir + os.sep + 'out.cpl')
        conf_lines = []
        l = f.readline()
        try:
            while l!='':
                conf_lines.append(
                    tuple(map(lambda x : float(x), l.rstrip().split(' ')))
                )
                l = f.readline()
        except ValueError:
            f.close()
            return []
        f.close()
        return conf_lines
        
    def run_pv(self, out_file='probs.dat', data_dir='.'):
        """Executes pv.

        out_file - Name of output file.
        data_dir - Where the data is found.
        """
        in_name = self._get_temp_file()
        out_name = self._get_temp_file()
        f = open(data_dir + os.sep + in_name, 'w')
        f.write('data_fst_outfile ' + out_file + ' out.dat\n')
        f.close()
        self._generate_intfile(data_dir)
        os.system('cd ' + data_dir + ' && ' +
                self._get_path('pv') + ' < ' + in_name + ' > ' + out_name)
        pvf = open(data_dir + os.sep + out_file, 'r')
        result = map(lambda x: tuple(map(lambda y: float(y), x.rstrip().split(' '))),
            pvf.readlines())
        pvf.close()
        os.remove(data_dir + os.sep + in_name)
        os.remove(data_dir + os.sep + out_name)
        return result

