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
import sys
from shutil import copyfile
from random import randint, random
from time import strftime, clock
#from logging import debug

if sys.version_info[0] == 3:
    maxint = sys.maxsize
else:
    maxint = sys.maxint

def my_float(f):
    #Because of Jython, mostly
    if f=="-nan": f="nan"
    return float(f)

class FDistController(object):
    def __init__(self, fdist_dir = '', ext = None):
        """Initializes the controller.
        
        fdist_dir is the directory where fdist2 is.
        ext is the extension of binaries (.exe on windows, 
          none on Unix)
        
        """
        self.tmp_idx = 0
        self.fdist_dir = fdist_dir
        self.os_name = os.name
        if sys.platform=='win32':
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

    def run_datacal(self, data_dir='.', version=1,
        crit_freq = 0.99, p = 0.5, beta= (0.25, 0.25)):
        """Executes datacal.
        
           data_dir - Where the data is found.
        """
        in_name = self._get_temp_file()
        out_name = self._get_temp_file()
        f = open(data_dir + os.sep + in_name, 'w')
        if version==1:
            f.write('a\n')
            datacal_name = "datacal"
        else:
            f.write('%f\n%f\n%f %f\na\n' % (crit_freq, p, beta[0], beta[1]))
            datacal_name = "Ddatacal"
        f.close()
        curr_dir = os.getcwd()
        os.system('cd ' + data_dir + ' && ' +
                self._get_path(datacal_name) + ' < ' + in_name + ' > ' + out_name)
        f = open(data_dir + os.sep + out_name)
        if version == 1:
            fst_line = f.readline().rstrip().split(' ')
            fst = my_float(fst_line[4])
            sample_line = f.readline().rstrip().split(' ')
            sample = int(sample_line[9])
        else:
            l = f.readline().rstrip().split(" ")
            loci, pops = int(l[-5]), int(l[-2])
            fst_line = f.readline().rstrip().split(' ')
            fst = my_float(fst_line[4])
            sample_line = f.readline().rstrip().split(' ')
            sample = int(sample_line[9])
            F_line = f.readline().rstrip().split(' ')
            F, obs = my_float(F_line[5]), int (F_line[8])
        f.close()
        os.remove(data_dir + os.sep + in_name)
        os.remove(data_dir + os.sep + out_name)
        if version==1:
            return fst, sample
        else:
            return fst, sample, loci, pops, F, obs

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
        mut = 0, num_sims = 50000, data_dir='.',
        is_dominant = False, theta = 0.06, beta = (0.25, 0.25),
        max_freq = 0.99):
        """Executes (d)fdist.
        
        Parameters:
        npops - Number of populations
        nsamples - Number of populations sampled
        fst - expected Fst
        sample_size - Sample size per population
                For dfdist: if zero a sample size file has to be provided 
        mut - 1=Stepwise, 0=Infinite allele
        num_sims - number of simulations
        data_dir - Where the data is found
        is_dominant - If true executes dfdist
        theta - Theta (=2Nmu)
        beta - Parameters for the beta prior
        max_freq - Maximum allowed frequency of the commonest allele

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
        if is_dominant:
            config_name = "Dfdist_params"
        else:
            config_name = "fdist_params2.dat"

        f = open(data_dir + os.sep + config_name, 'w')
        f.write(str(npops) + '\n')
        f.write(str(nsamples) + '\n')
        f.write(str(fst) + '\n')
        f.write(str(sample_size) + '\n')
        if is_dominant:
            f.write(str(theta) + '\n')
        else:
            f.write(str(mut) + '\n')
        f.write(str(num_sims) + '\n')
        if is_dominant:
            f.write("%f %f\n" % beta)
            f.write("%f\n" % max_freq)
        f.close()
        self._generate_intfile(data_dir)

        if is_dominant:
            bin_name = "Dfdist"
        else:
            bin_name = "fdist2"
        os.system('cd ' + data_dir + ' && ' +
            self._get_path(bin_name) + ' < ' + in_name + ' > ' + out_name)
        f = open(data_dir + os.sep + out_name)
        lines = f.readlines()
        f.close()
        for line in lines:
          if line.startswith('average Fst'):
            fst = my_float(line.rstrip().split(' ')[-1])
        os.remove(data_dir + os.sep + in_name)
        os.remove(data_dir + os.sep + out_name)
        return fst

    def run_fdist_force_fst(self, npops, nsamples, fst, sample_size,
        mut = 0, num_sims = 50000, data_dir='.',
        try_runs = 5000, limit=0.001,
        is_dominant = False, theta = 0.06, beta = (0.25, 0.25),
        max_freq = 0.99):
        """Executes fdist trying to force Fst.
        
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
            real_fst = self.run_fdist(npops, nsamples,
                current_run_fst, sample_size,
                mut, try_runs, data_dir,
                is_dominant, theta, beta, max_freq)
            #debug('got real fst ' +  str(real_fst))
            if abs(real_fst - fst) < limit:
                #debug('We are OK')
                return self.run_fdist(npops, nsamples, current_run_fst,
                    sample_size,
                    mut, num_sims, data_dir,
                    is_dominant, theta, beta, max_freq)
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
                        sample_size, mut, num_sims, data_dir,
                        is_dominant, theta, beta, max_freq)
                current_run_fst = (max_run_fst + current_run_fst)/2

    def run_cplot(self, ci= 0.95, data_dir='.', version = 1, smooth=0.04):
        """Executes cplot.

        ci - Confidence interval.
        data_dir - Where the data is found.
        """
        in_name = self._get_temp_file()
        out_name = self._get_temp_file()
        f = open(data_dir + os.sep + in_name, 'w')
        if version == 1:
            f.write('out.dat out.cpl\n' + str(ci) + '\n')
        else:
            f.write("\n".join([
                "data_fst_outfile out.cpl out.dat",
                str(ci), str(smooth)]))
        f.close()
        curr_dir = os.getcwd()
        self._generate_intfile(data_dir)
        if version == 1:
            cplot_name = "cplot"
        else:
            cplot_name = "cplot2"
        os.system('cd ' + data_dir + ' && '  +
            self._get_path(cplot_name) + ' < ' + in_name + ' > ' + out_name)
        os.remove(data_dir + os.sep + in_name)
        os.remove(data_dir + os.sep + out_name)
        f = open(data_dir + os.sep + 'out.cpl')
        conf_lines = []
        l = f.readline()
        try:
            while l!='':
                conf_lines.append(
                    tuple(map(lambda x : my_float(x), l.rstrip().split(' ')))
                )
                l = f.readline()
        except ValueError:
            f.close()
            return []
        f.close()
        return conf_lines
        
    def run_pv(self, out_file='probs.dat', data_dir='.',
               version = 1, smooth=0.04):
        """Executes pv.

        out_file - Name of output file.
        data_dir - Where the data is found.
        """
        in_name = self._get_temp_file()
        out_name = self._get_temp_file()
        f = open(data_dir + os.sep + in_name, 'w')
        f.write('data_fst_outfile ' + out_file + ' out.dat\n')
        f.write(str(smooth) + '\n')
        f.close()
        self._generate_intfile(data_dir)
        if version == 1:
            pv_name = "pv"
        else:
            pv_name = "pv2"
        os.system('cd ' + data_dir + ' && ' +
                self._get_path(pv_name) + ' < ' + in_name + ' > ' + out_name)
        pvf = open(data_dir + os.sep + out_file, 'r')
        result = map(lambda x: tuple(map(lambda y: my_float(y), x.rstrip().split(' '))),
            pvf.readlines())
        pvf.close()
        os.remove(data_dir + os.sep + in_name)
        os.remove(data_dir + os.sep + out_name)
        return result

