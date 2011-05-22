# Copyright 2009 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.



"""
This module allows to control GenePop.

"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

from Bio.Application import AbstractCommandline, _Argument, _Option

def _gp_float(tok):
    """Gets a float from a token, if it fails, returns the string.
    """
    try:
        return float(tok)
    except ValueError:
        return str(tok)

def _gp_int(tok):
    """Gets a int from a token, if it fails, returns the string.
    """
    try:
        return int(tok)
    except ValueError:
        return str(tok)
    
    
def _read_allele_freq_table(f):
    l = f.readline()
    while l.find(" --")==-1:
        if l == "":
            raise StopIteration
        if l.find("No data")>-1:
            return None, None
        l = f.readline()
    alleles = filter(lambda x: x != '', f.readline().rstrip().split(" "))
    alleles = map(lambda x: _gp_int(x), alleles)
    l = f.readline().rstrip()
    table = []
    while l != "":
        line = filter(lambda x: x != '', l.split(" "))
        try:
            table.append(
                (line[0],
                map(lambda x: _gp_float(x), line[1:-1]),
                _gp_int(line[-1])))
        except ValueError:
            table.append(
                (line[0],
                [None] * len(alleles),
                0))
        l = f.readline().rstrip()
    return alleles, table

def _read_table(f, funs):
    table = []
    l = f.readline().rstrip()
    while l.find("---")==-1:
        l = f.readline().rstrip()
    l = f.readline().rstrip()
    while l.find("===")==-1 and l.find("---")==-1 and l != "":
        toks = filter(lambda x: x != "", l.split(" ")) 
        line = []
        for i in range(len(toks)):
            try:
                line.append(funs[i](toks[i]))
            except ValueError:
                line.append(toks[i])  # Could not cast
        table.append(tuple(line))
        l = f.readline().rstrip()
    return table

def _read_triangle_matrix(f):
    matrix = []
    l = f.readline().rstrip()
    while l != "":
        matrix.append(
            map(lambda x: _gp_float(x),
                filter(lambda y: y != "", l.split(" "))))
        l = f.readline().rstrip()
    return matrix

def _read_headed_triangle_matrix(f):
    matrix = {}
    header = f.readline().rstrip()
    if header.find("---")>-1 or header.find("===")>-1:
        header = f.readline().rstrip()
    nlines = len(filter(lambda x:x != '', header.split(' '))) - 1
    for line_pop in range(nlines):
        l = f.readline().rstrip()
        vals = filter(lambda x:x != '', l.split(' ')[1:])
        clean_vals = []
        for val in vals:
            try:
                clean_vals.append(_gp_float(val))
            except ValueError:
                clean_vals.append(None)
        for col_pop in range(len(clean_vals)):
            matrix[(line_pop+1, col_pop)] = clean_vals[col_pop]
    return matrix

def _hw_func(stream, is_locus, has_fisher = False): 
    l = stream.readline()
    if is_locus:
        hook = "Locus "
    else:
        hook = " Pop : "
    while l != "":
        if l.startswith(hook):
            stream.readline()
            stream.readline()
            stream.readline()
            table = _read_table(stream,[str,_gp_float,_gp_float,_gp_float,_gp_float,_gp_int,str])
            #loci might mean pop if hook="Locus "
            loci = {}
            for entry in table:
                if len(entry) < 3:
                    loci[entry[0]] = None
                else:
                    locus, p, se, fis_wc, fis_rh, steps = entry[:-1]
                    if se == "-": se = None
                    loci[locus] = p, se, fis_wc, fis_rh, steps
            return loci
        l = stream.readline()
    #self.done = True
    raise StopIteration

class _FileIterator:
    """Iterator which crawls over a stream of lines with a function.

       The generator function is expected to yield a tuple, while
       consuming input
    """
    def __init__(self, func, stream, fname):
        self.func = func
        self.stream = stream
        self.fname = fname
        self.done = False

    def __iter__(self):
        if self.done:
            self.done = True
            raise StopIteration
        return self

    def next(self):
        return self.func(self)

    def __del__(self):
        self.stream.close()
        try:
            os.remove(self.fname)
        except OSError:
            #Jython seems to call the iterator twice
            pass

class _GenePopCommandline(AbstractCommandline):
    """ Command Line Wrapper for GenePop.
    """
    def __init__(self, genepop_dir=None, cmd='Genepop', **kwargs):
        self.parameters = [
                _Argument(["command"],
                    "GenePop option to be called",
                    is_required=True),
                _Argument(["mode"],
                    "Should allways be batch",
                    is_required=True),
                _Argument(["input"],
                    "Input file",
                    is_required=True),
                _Argument(["Dememorization"],
                    "Dememorization step"),
                _Argument(["BatchNumber"],
                    "Number of MCMC batches"),
                _Argument(["BatchLength"],
                    "Length of MCMC chains"),
                _Argument(["HWtests"],
                    "Enumeration or MCMC"),
                _Argument(["IsolBDstatistic"],
                    "IBD statistic (a or e)"),
                _Argument(["MinimalDistance"],
                    "Minimal IBD distance"),
                _Argument(["GeographicScale"],
                    "Log or Linear"),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
        self.set_parameter("mode", "Mode=Batch")

    def set_menu(self, option_list):
        """Sets the menu option.

        Example set_menu([6,1]) = get all F statistics (menu 6.1)
        """
        self.set_parameter("command", "MenuOptions="+
                ".".join(map(lambda x:str(x),option_list)))

    def set_input(self, fname):
        """Sets the input file name.
        """
        self.set_parameter("input", "InputFile="+fname)

class GenePopController(object):
    def __init__(self, genepop_dir = None):
        """Initializes the controller.
        
        genepop_dir is the directory where GenePop is.

        The binary should be called Genepop (capital G)
        
        """
        self.controller = _GenePopCommandline(genepop_dir)

    def _remove_garbage(self, fname_out):
        try:
            if fname_out != None: os.remove(fname_out)
        except OSError:
            pass # safe
        try:
            os.remove("genepop.txt")
        except OSError:
            pass # safe
        try:
            os.remove("fichier.in")
        except OSError:
            pass # safe
        try:
            os.remove("cmdline.txt")
        except OSError:
            pass # safe

    def _get_opts(self, dememorization, batches, iterations, enum_test=None):
        opts = {}
        opts["Dememorization"]=dememorization
        opts["BatchNumber"]=batches
        opts["BatchLength"]=iterations
        if enum_test != None:
            if enum_test == True:
                opts["HWtests"]="Enumeration"
            else:
                opts["HWtests"]="MCMC"
        return opts

    def _run_genepop(self, extensions, option, fname, opts={}):
        for extension in extensions:
            self._remove_garbage(fname + extension)
        self.controller.set_menu(option)
        self.controller.set_input(fname)
        for opt in opts:
            self.controller.set_parameter(opt, opt+"="+str(opts[opt]))
        self.controller() #checks error level is zero
        self._remove_garbage(None)
        return


    def _test_pop_hz_both(self, fname, type, ext, enum_test = True,
        dememorization = 10000, batches = 20, iterations = 5000):
        """Hardy-Weinberg test for heterozygote deficiency/excess.

           Returns a population iterator containg
               A dictionary[locus]=(P-val, SE, Fis-WC, Fis-RH, steps)
                 Some loci have a None if the info is not available
                 SE might be none (for enumerations)
        """
        opts = self._get_opts(dememorization, batches, iterations, enum_test)
        self._run_genepop([ext], [1, type], fname, opts)
        f = open(fname + ext)
        def hw_func(self):
            return _hw_func(self.stream, False)
        return _FileIterator(hw_func, f, fname + ext)

    def _test_global_hz_both(self, fname, type, ext, enum_test = True,
        dememorization = 10000, batches = 20, iterations = 5000):
        """Global Hardy-Weinberg test for heterozygote deficiency/excess.

           Returns a triple with:
             A list per population containg
               (pop_name, P-val, SE, switches)
                 Some pops have a None if the info is not available
                 SE might be none (for enumerations)
             A list per loci containg
               (locus_name, P-val, SE, switches)
                 Some loci have a None if the info is not available
                 SE might be none (for enumerations)
             Overall results (P-val, SE, switches)

        """
        opts = self._get_opts(dememorization, batches, iterations, enum_test)
        self._run_genepop([ext], [1, type], fname, opts)
        def hw_pop_func(self):
            return _read_table(self.stream, [str, _gp_float, _gp_float, _gp_float])
        f1 = open(fname + ext)
        l = f1.readline()
        while l.find("by population") == -1:
            l = f1.readline()
        pop_p = _read_table(f1, [str, _gp_float, _gp_float, _gp_float])
        f2 = open(fname + ext)
        l = f2.readline()
        while l.find("by locus") == -1:
            l = f2.readline()
        loc_p = _read_table(f2, [str, _gp_float, _gp_float, _gp_float])
        f = open(fname + ext)
        l = f.readline()
        while l.find("all locus") == -1:
            l = f.readline()
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        l = f.readline().rstrip()
        p, se, switches = tuple(map(lambda x: _gp_float(x),
            filter(lambda y: y != "",l.split(" "))))
        f.close()
        return pop_p, loc_p, (p, se, switches)

    #1.1
    def test_pop_hz_deficiency(self, fname, enum_test = True,
        dememorization = 10000, batches = 20, iterations = 5000):
        """Hardy-Weinberg test for heterozygote deficiency.

           Returns a population iterator containg
               A dictionary[locus]=(P-val, SE, Fis-WC, Fis-RH, steps)
                 Some loci have a None if the info is not available
                 SE might be none (for enumerations)
        """
        return self._test_pop_hz_both(fname, 1, ".D", enum_test,
            dememorization, batches, iterations)

    #1.2
    def test_pop_hz_excess(self, fname, enum_test = True,
        dememorization = 10000, batches = 20, iterations = 5000):
        """Hardy-Weinberg test for heterozygote deficiency.

           Returns a population iterator containg
               A dictionary[locus]=(P-val, SE, Fis-WC, Fis-RH, steps)
                 Some loci have a None if the info is not available
                 SE might be none (for enumerations)
        """
        return self._test_pop_hz_both(fname, 2, ".E", enum_test,
            dememorization, batches, iterations)

    #1.3 P file
    def test_pop_hz_prob(self, fname, ext, enum_test = False,
        dememorization = 10000, batches = 20, iterations = 5000):
        """Hardy-Weinberg test based on probability.

           Returns 2 iterators and a final tuple:

          1. Returns a loci iterator containing
               b. A dictionary[pop_pos]=(P-val, SE, Fis-WC, Fis-RH, steps)
                 Some pops have a None if the info is not available
                 SE might be none (for enumerations)
               c. Result of Fisher's test (Chi2, deg freedom, prob)
          2. Returns a population iterator containg
               a. A dictionary[locus]=(P-val, SE, Fis-WC, Fis-RH, steps)
                 Some loci have a None if the info is not available
                 SE might be none (for enumerations)
               b. Result of Fisher's test (Chi2, deg freedom, prob)
          3. (Chi2, deg freedom, prob)
        """
        opts = self._get_opts(dememorization, batches, iterations, enum_test)
        self._run_genepop([ext], [1, 3], fname, opts)
        def hw_prob_loci_func(self): 
            return  _hw_func(self.stream, True, True)
        def hw_prob_pop_func(self): 
            return _hw_func(self.stream, False, True)
        shutil.copyfile(fname+".P", fname+".P2")
        f1 = open(fname + ".P")
        f2 = open(fname + ".P2")
        return _FileIterator(hw_prob_loci_func, f1, fname + ".P"), _FileIterator(hw_prob_pop_func, f2, fname + ".P2")

    #1.4
    def test_global_hz_deficiency(self, fname, enum_test = True,
        dememorization = 10000, batches = 20, iterations = 5000):
        """Global Hardy-Weinberg test for heterozygote deficiency.

           Returns a triple with:
             An list per population containg
               (pop_name, P-val, SE, switches)
                 Some pops have a None if the info is not available
                 SE might be none (for enumerations)
             An list per loci containg
               (locus_name, P-val, SE, switches)
                 Some loci have a None if the info is not available
                 SE might be none (for enumerations)
             Overall results (P-val, SE, switches)
        """
        return self._test_global_hz_both(fname, 4, ".DG", enum_test,
            dememorization, batches, iterations)


    #1.5
    def test_global_hz_excess(self, fname, enum_test = True,
        dememorization = 10000, batches = 20, iterations = 5000):
        """Global Hardy-Weinberg test for heterozygote excess.

           Returns a triple with:
             An list per population containg
               (pop_name, P-val, SE, switches)
                 Some pops have a None if the info is not available
                 SE might be none (for enumerations)
             An list per loci containg
               (locus_name, P-val, SE, switches)
                 Some loci have a None if the info is not available
                 SE might be none (for enumerations)
             Overall results (P-val, SE, switches)
        """
        return self._test_global_hz_both(fname, 5, ".EG", enum_test,
            dememorization, batches, iterations)

    #2.1
    def test_ld(self, fname,
        dememorization = 10000, batches = 20, iterations = 5000):
        opts = self._get_opts(dememorization, batches, iterations)
        self._run_genepop([".DIS"], [2, 1], fname, opts)
        def ld_pop_func(self):
            current_pop = None
            l = self.stream.readline().rstrip()
            if l == "":
                self.done = True
                raise StopIteration
            toks = filter(lambda x: x != "", l.split(" "))
            pop, locus1, locus2 = toks[0], toks[1], toks[2]
            if not hasattr(self, "start_locus1"):
                start_locus1, start_locus2 = locus1, locus2
                current_pop = -1
            if locus1 == start_locus1 and locus2 == start_locus2:
                current_pop += 1
            if toks[3] == "No":
                return current_pop, pop, (locus1, locus2), None
            p, se, switches = _gp_float(toks[3]), _gp_float(toks[4]), _gp_int(toks[5])
            return current_pop, pop, (locus1, locus2), (p, se, switches)
        def ld_func(self):
            l = self.stream.readline().rstrip()
            if l == "":
                self.done = True
                raise StopIteration
            toks = filter(lambda x: x != "", l.split(" "))
            locus1, locus2 = toks[0], toks[2]
            try:
                chi2, df, p = _gp_float(toks[3]), _gp_int(toks[4]), _gp_float(toks[5])
            except ValueError:
                return (locus1, locus2), None
            return (locus1, locus2), (chi2, df, p)
        f1 = open(fname + ".DIS")
        l = f1.readline()
        while l.find("----")==-1:
            l = f1.readline()
        shutil.copyfile(fname + ".DIS", fname + ".DI2")
        f2 = open(fname + ".DI2")
        l = f2.readline()
        while l.find("Locus pair")==-1:
            l = f2.readline()
        while l.find("----")==-1:
            l = f2.readline()
        return _FileIterator(ld_pop_func, f1, fname+".DIS"), _FileIterator(ld_func, f2, fname + ".DI2")

    #2.2
    def create_contingency_tables(self, fname):
        raise NotImplementedError

    #3.1 PR/GE files
    def test_genic_diff_all(self, fname,
        dememorization = 10000, batches = 20, iterations = 5000):
        raise NotImplementedError

    #3.2 PR2/GE2 files
    def test_genic_diff_pair(self, fname,
        dememorization = 10000, batches = 20, iterations = 5000):
        raise NotImplementedError

    #3.3 G files
    def test_genotypic_diff_all(self, fname,
        dememorization = 10000, batches = 20, iterations = 5000):
        raise NotImplementedError

    #3.4 2G2 files
    def test_genotypic_diff_pair(self, fname,
        dememorization = 10000, batches = 20, iterations = 5000):
        raise NotImplementedError

    #4
    def estimate_nm(self, fname):
        self._run_genepop(["PRI"], [4], fname)
        f = open(fname + ".PRI")
        lines = f.readlines() # Small file, it is ok
        f.close()
        for line in lines:
            m = re.search("Mean sample size: ([.0-9]+)", line)
            if m != None:
                mean_sample_size = _gp_float(m.group(1))
            m = re.search("Mean frequency of private alleles p\(1\)= ([.0-9]+)", line)
            if m != None:
                mean_priv_alleles = _gp_float(m.group(1))
            m = re.search("N=10: ([.0-9]+)", line)
            if m != None:
                mig10 = _gp_float(m.group(1))
            m = re.search("N=25: ([.0-9]+)", line)
            if m != None:
                 mig25 = _gp_float(m.group(1))
            m = re.search("N=50: ([.0-9]+)", line)
            if m != None:
                 mig50 = _gp_float(m.group(1))
            m = re.search("for size= ([.0-9]+)", line)
            if m != None:
                 mig_corrected = _gp_float(m.group(1))
        os.remove(fname + ".PRI")
        return mean_sample_size, mean_priv_alleles, mig10, mig25, mig50, mig_corrected

    #5.1
    def calc_allele_genotype_freqs(self, fname):
        """Calculates allele and genotype frequencies per locus and per sample.
        
        Parameters:
        fname - file name

        Returns tuple with 2 elements:
        Population iterator with
            population name
            Locus dictionary with key = locus name and content tuple as
              Genotype List with
                (Allele1, Allele2, observed, expected)
              (expected homozygotes, observed hm, 
              expected heterozygotes, observed ht)
              Allele frequency/Fis dictionary with allele as key and
                (count, frequency, Fis Weir & Cockerham)
              Totals as a pair
                count
                Fis Weir & Cockerham,
                Fis Robertson & Hill
        Locus iterator with
            Locus name
            allele list
            Population list with a triple
               population name
               list of allele frequencies in the same order as allele list above
               number of genes

        
        Will create a file called fname.INF
        """
        self._run_genepop(["INF"], [5,1], fname)
        #First pass, general information
        #num_loci = None
        #num_pops = None
        #f = open(fname + ".INF")
        #l = f.readline()
        #while (num_loci == None or num_pops == None) and l != '':
        #    m = re.search("Number of populations detected : ([0-9+])", l)
        #    if m != None:
        #        num_pops = _gp_int(m.group(1))
        #    m = re.search("Number of loci detected        : ([0-9+])", l)
        #    if m != None:
        #        num_loci = _gp_int(m.group(1))
        #    l = f.readline()
        #f.close()
        def pop_parser(self):
            if hasattr(self, "old_line"):
                l = self.old_line
                del self.old_line
            else:
                l = self.stream.readline()
            loci_content = {}
            while l != '':
                l = l.rstrip()
                if l.find("Tables of allelic frequencies for each locus")>-1:
                    return self.curr_pop, loci_content
                match = re.match(".*Pop: (.+) Locus: (.+)", l)
                if match != None:
                    pop = match.group(1)
                    locus = match.group(2)
                    if not hasattr(self, "first_locus"):
                        self.first_locus = locus
                    if hasattr(self, "curr_pop"):
                        if self.first_locus == locus:
                            old_pop = self.curr_pop
                            #self.curr_pop = pop
                            self.old_line = l
                            del self.first_locus
                            del self.curr_pop
                            return old_pop, loci_content
                    self.curr_pop = pop
                else:
                    l = self.stream.readline()
                    continue
                geno_list = []
                l = self.stream.readline()
                if l.find("No data")>-1: continue

                while l.find("Genotypes  Obs.")==-1:
                    l = self.stream.readline()

                while l != "\n":
                    m2 = re.match(" +([0-9]+) , ([0-9]+) *([0-9]+) *(.+)",l)
                    if m2 != None:
                        geno_list.append((_gp_int(m2.group(1)), _gp_int(m2.group(2)),
                            _gp_int(m2.group(3)), _gp_float(m2.group(4))))
                    else:
                        l = self.stream.readline()
                        continue
                    l = self.stream.readline()

                while l.find("Expected number of ho")==-1:
                    l = self.stream.readline()
                expHo =  _gp_float(l[38:])
                l = self.stream.readline()
                obsHo =  _gp_int(l[38:])
                l = self.stream.readline()
                expHe =  _gp_float(l[38:])
                l = self.stream.readline()
                obsHe =  _gp_int(l[38:])
                l = self.stream.readline()

                while l.find("Sample count")==-1:
                    l = self.stream.readline()
                l = self.stream.readline()
                freq_fis={}
                overall_fis = None
                while l.find("----")==-1:
                    vals = filter(lambda x: x!='',
                            l.rstrip().split(' '))
                    if vals[0]=="Tot":
                        overall_fis = _gp_int(vals[1]), \
                                _gp_float(vals[2]), _gp_float(vals[3])
                    else:
                        freq_fis[_gp_int(vals[0])] = _gp_int(vals[1]), \
                                _gp_float(vals[2]), _gp_float(vals[3])
                    l = self.stream.readline()
                loci_content[locus] = geno_list, \
                        (expHo, obsHo, expHe, obsHe), \
                        freq_fis, overall_fis
            self.done = True
            raise StopIteration
        def locus_parser(self):
            l = self.stream.readline()
            while l != "":
                l = l.rstrip()
                match = re.match(" Locus: (.+)", l)
                if match != None:
                    locus = match.group(1)
                    alleles, table = _read_allele_freq_table(self.stream)
                    return locus, alleles, table
                l = self.stream.readline()
            self.done = True
            raise StopIteration

        popf = open(fname + ".INF")
        shutil.copyfile(fname + ".INF", fname + ".IN2")
        locf = open(fname + ".IN2")
        pop_iter = _FileIterator(pop_parser, popf, fname + ".INF")
        locus_iter = _FileIterator(locus_parser, locf, fname + ".IN2")
        return (pop_iter, locus_iter)

    def _calc_diversities_fis(self, fname, ext):
        self._run_genepop([ext], [5,2], fname)
        f = open(fname + ext)
        l = f.readline()
        while l != "":
            l = l.rstrip()
            if l.startswith("Statistics per sample over all loci with at least two individuals typed"):
                avg_fis = _read_table(f, [str, _gp_float, _gp_float, _gp_float])
                avg_Qintra = _read_table(f, [str, _gp_float])
            l = f.readline()
        f.close()
        def fis_func(self): 
            l = self.stream.readline()
            while l != "":
                l = l.rstrip()
                m = re.search("Locus: (.+)", l)
                if m != None:
                    locus = m.group(1)
                    self.stream.readline()
                    if self.stream.readline().find("No complete")>-1: return locus, None
                    self.stream.readline()
                    fis_table = _read_table(self.stream, [str, _gp_float, _gp_float, _gp_float])
                    self.stream.readline()
                    avg_qinter, avg_fis = tuple(map (lambda x: _gp_float(x),
                        filter(lambda y:y != "", self.stream.readline().split(" "))))
                    return locus, fis_table, avg_qinter, avg_fis
                l = self.stream.readline()
            self.done = True
            raise StopIteration
        dvf = open(fname + ext)
        return _FileIterator(fis_func, dvf, fname + ext), avg_fis, avg_Qintra

    #5.2
    def calc_diversities_fis_with_identity(self, fname):
        return self._calc_diversities_fis(fname, ".DIV")

    #5.3
    def calc_diversities_fis_with_size(self, fname):
        raise NotImplementedError

    #6.1 Less genotype frequencies
    def calc_fst_all(self, fname):
        """Executes GenePop and gets Fst/Fis/Fit (all populations)
        
        Parameters:
        fname - file name

        Returns:
        (multiLocusFis, multiLocusFst, multiLocus Fit),
        Iterator of tuples
          (Locus name, Fis, Fst, Fit, Qintra, Qinter)

        Will create a file called fname.FST .

        This does not return the genotype frequencies.
        
        """
        self._run_genepop([".FST"], [6,1], fname)
        f = open(fname + ".FST")
        l = f.readline()
        while l != '':
            if l.startswith('           All:'):
                toks=filter(lambda x:x!="", l.rstrip().split(' '))
                try:
                    allFis = _gp_float(toks[1])
                except ValueError:
                    allFis = None
                try:
                    allFst = _gp_float(toks[2])
                except ValueError:
                    allFst = None
                try:
                    allFit = _gp_float(toks[3])
                except ValueError:
                    allFit = None
            l = f.readline()
        f.close()
        f = open(fname + ".FST")
        def proc(self): 
            if hasattr(self, "last_line"):
                l = self.last_line
                del self.last_line
            else:
                l = self.stream.readline()
            locus = None
            fis = None
            fst = None
            fit = None
            qintra = None
            qinter = None
            while l != '':
                l = l.rstrip()
                if l.startswith('  Locus:'):
                    if locus != None:
                        self.last_line = l
                        return locus, fis, fst, fit, qintra, qinter
                    else:
                        locus = l.split(':')[1].lstrip()
                elif l.startswith('Fis^='):
                    fis = _gp_float(l.split(' ')[1])
                elif l.startswith('Fst^='):
                    fst = _gp_float(l.split(' ')[1])
                elif l.startswith('Fit^='):
                    fit = _gp_float(l.split(' ')[1])
                elif l.startswith('1-Qintra^='):
                    qintra = _gp_float(l.split(' ')[1])
                elif l.startswith('1-Qinter^='):
                    qinter = _gp_float(l.split(' ')[1])
                    return locus, fis, fst, fit, qintra, qinter
                l = self.stream.readline()
            if locus != None:
                return locus, fis, fst, fit, qintra, qinter
            self.stream.close()
            self.done = True
            raise StopIteration
        return (allFis, allFst, allFit), _FileIterator(proc , f, fname + ".FST")

    #6.2
    def calc_fst_pair(self, fname):
        self._run_genepop([".ST2", ".MIG"], [6,2], fname)
        f = open(fname + ".ST2")
        l = f.readline()
        while l != "":
            l = l.rstrip()
            if l.startswith("Estimates for all loci"):
                avg_fst = _read_headed_triangle_matrix(f)
            l = f.readline()
        f.close()
        def loci_func(self): 
            l = self.stream.readline()
            while l != "":
                l = l.rstrip()
                m = re.search(" Locus: (.+)", l)
                if m != None:
                    locus = m.group(1)
                    matrix = _read_headed_triangle_matrix(self.stream)
                    return locus, matrix
                l = self.stream.readline()
            self.done = True
            raise StopIteration
        stf = open(fname + ".ST2")
        os.remove(fname + ".MIG")
        return _FileIterator(loci_func, stf, fname + ".ST2"), avg_fst

    #6.3
    def calc_rho_all(self, fname):
        raise NotImplementedError

    #6.4
    def calc_rho_pair(self, fname):
        raise NotImplementedError

    def _calc_ibd(self, fname, sub, stat="a", scale="Log", min_dist=0.00001):
        """Calculates isolation by distance statistics
        """
        self._run_genepop([".GRA", ".MIG", ".ISO"], [6,sub],
            fname, opts = {
            "MinimalDistance" : min_dist,
            "GeographicScale" : scale,
            "IsolBDstatistic" : stat,
            })
        f = open(fname + ".ISO")
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        estimate = _read_triangle_matrix(f)
        f.readline()
        f.readline()
        distance = _read_triangle_matrix(f)
        f.readline()
        match = re.match("a = (.+), b = (.+)", f.readline().rstrip())
        a = _gp_float(match.group(1))
        b = _gp_float(match.group(2))
        f.readline()
        f.readline()
        match = re.match(" b=(.+)", f.readline().rstrip())
        bb = _gp_float(match.group(1))
        match = re.match(".*\[(.+)  ;  (.+)\]", f.readline().rstrip())
        bblow = _gp_float(match.group(1))
        bbhigh = _gp_float(match.group(2))
        f.close()
        os.remove(fname + ".MIG")
        os.remove(fname + ".GRA")
        os.remove(fname + ".ISO")
        return estimate, distance, (a, b), (bb, bblow, bbhigh)

    #6.5
    def calc_ibd_diplo(self, fname, stat="a", scale="Log", min_dist=0.00001):
        """Calculates isolation by distance statistics for diploid data.

           See _calc_ibd for parameter details.
           Note that each pop can only have a single individual and
           the individual name has to be the sample coordinates.
        """
        return self._calc_ibd(fname, 5, stat, scale, min_dist)

    #6.6
    def calc_ibd_haplo(self, fname, stat="a", scale="Log", min_dist=0.00001):
        """Calculates isolation by distance statistics for haploid data.

           See _calc_ibd for parameter details.
           Note that each pop can only have a single individual and
           the individual name has to be the sample coordinates.
        """
        return self._calc_ibd(fname, 6, stat, scale, min_dist)


