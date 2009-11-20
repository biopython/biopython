# Copyright 2007 by Tiago Antao.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from math import sqrt
from sys import argv,exit
from os import sep, mkdir
import re

from Bio.PopGen.SimCoal import builtin_tpl_dir


def exec_template(template):
    executed_template = template
    match = re.search('!!!(.*?)!!!', executed_template, re.MULTILINE)
    #while len(match.groups())>0:
    while match:
        exec_result = str(eval(match.groups()[0]))
        executed_template = executed_template.replace(
               '!!!' + match.groups()[0] + '!!!',
               exec_result, 1)
        match = re.search('!!!(.*?)!!!', executed_template, re.MULTILINE)
        #match = patt.matcher(String(executed_template))
    return executed_template
    

def process_para(in_string, out_file_prefix, para_list, curr_values):
    if (para_list == []):
        template = in_string
        f_name = out_file_prefix
        #f_name += '_' + str(total_size)
        for tup in curr_values:
            name, val = tup
            f_name += '_' + str(val) 
            #reg = re.compile('\?' + name, re.MULTILINE)
            #template = re.sub(reg, str(val), template)
            template = template.replace('?'+name, str(val))
        f = open(f_name + '.par', 'w')
        #executed_template = template
        executed_template = exec_template(template)
        clean_template =  executed_template.replace('\r\n','\n').replace('\n\n','\n')
        f.write(clean_template)
        f.close()
        return [f_name]
    else:
        name, rng = para_list[0]
        fnames = []
        for val in rng:
            new_values = [(name, val)]
            new_values.extend(curr_values)
            more_names = process_para(in_string, out_file_prefix, para_list[1:], new_values)
            fnames.extend(more_names)
        return fnames


def dupe(motif, times):
    ret_str = ''
    for i in range(1, times + 1):
        ret_str += motif + '\r\n'
    return ret_str

def get_xy_from_matrix(x_max, y_max, pos):
    y = (pos-1) / x_max
    x = (pos-1) % x_max
    return x, y

def get_step_2d(x_max, y_max, x, y, mig):
    my_x,    my_y    = get_xy_from_matrix(x_max, y_max, y)
    other_x, other_y = get_xy_from_matrix(x_max, y_max, x)
        
    if (my_x-other_x)**2 + (my_y-other_y)**2 == 1:
        return str(mig) + ' '
    else:
        return '0 '

def generate_ssm2d_mat(x_max, y_max, mig):
    mig_mat = ''
    for x in range(1, x_max*y_max + 1):
        for y in range(1, x_max*y_max + 1):
            mig_mat += get_step_2d(x_max, y_max, x, y, mig)
        mig_mat += "\r\n"
    return mig_mat

def generate_island_mat(total_size, mig):
    mig_mat = ''
    for x in range(1, total_size + 1):
        for y in range(1, total_size + 1):
            if (x == y):
                mig_mat += '0 '
            else:
                mig_mat += '!!!' + str(mig) + '!!! '
        mig_mat += "\r\n"
    return mig_mat

def generate_null_mat(total_size):
    null_mat = ''
    for x in range(1, total_size + 1):
        for y in range(1, total_size + 1):
            null_mat += '0 '
        null_mat += '\r\n'
    return null_mat

def generate_join_events(t, total_size, join_size, orig_size):
    events = ''
    for i in range(1, total_size-1):
        events += str(t) + ' ' + str(i) + ' 0 1 1 0 1\r\n'
    events += str(t) + ' ' + str(total_size-1) + ' 0 1 ' + str(1.0*total_size*join_size/orig_size) + ' 0 1\r\n'
    return events

def no_processor(in_string):
    return in_string

def process_text(in_string, out_file_prefix, para_list, curr_values,
                 specific_processor):
    text = specific_processor(in_string)
    return process_para(text, out_file_prefix, para_list, [])

#def prepare_dir():
#    try:
#        mkdir(sep.join([Config.dataDir, 'SimCoal'])) #Should exist, but...
#    except OSError:
#        pass #Its ok if already exists
#    try:
#        mkdir(sep.join([Config.dataDir, 'SimCoal', 'runs']))
#    except OSError:
#        pass #Its ok if already exists
    

#sep is because of jython
def generate_model(par_stream, out_prefix, params,
    specific_processor = no_processor, out_dir = '.'):
    #prepare_dir()
    text = par_stream.read()
    out_file_prefix = sep.join([out_dir, out_prefix])
    return process_text(text, out_file_prefix, params, [], specific_processor)


def get_demography_template(stream, model, tp_dir = None):
    '''
        Gets a demograpy template.
 
        Most probably this model needs to be sent to GenCases.
 
        stream - Writable stream.
        param  - Template file.
        tp_dir - Directory where to find the template, if None
                 use an internal template
    '''
    if tp_dir == None:
        #Internal Template
        f = open(sep.join([builtin_tpl_dir, model + '.par']), 'r')
    else:
        #External template
        f = open(sep.join([tp_dir, model + '.par']), 'r')
    l = f.readline()
    while l!='':
        stream.write(l)
        l = f.readline()
    f.close()

def _gen_loci(stream, loci):
    stream.write('//Number of contiguous linkage blocks in chromosome\n')
    stream.write(str(len(loci)) + '\n')
    stream.write('//Per Block: Data type, No. of loci, Recombination rate to the right-side locus, plus optional parameters\n')
    for locus in loci:
        stream.write(' '.join([locus[0]] +
            map(lambda x: str(x), list(locus[1])
        )) + '\n')

def get_chr_template(stream, chrs):
    '''
        Writes a Simcoal2 loci template part.

        stream - Writable stream.
        chr    - Chromosome list.

        Current loci list:
          [(chr_repeats,[(marker, (params))])]
          chr_repeats --> Number of chromosome repeats
          marker  --> 'SNP', 'DNA', 'RFLP', 'MICROSAT'
          params  --> Simcoal2 parameters for markers (list of floats
            or ints - if to be processed by generate_model)
    '''
    num_chrs = reduce(lambda x, y: x + y[0], chrs, 0)
    stream.write('//Number of independent (unlinked) chromosomes, and "chromosome structure" flag:  0 for identical structure across chromosomes, and  1 for different structures on different chromosomes.\n')
    if len(chrs) > 1 or num_chrs == 1:
        stream.write(str(num_chrs) + ' 1\n')
    else:
        stream.write(str(num_chrs) + ' 0\n')
    for chr in chrs:
        repeats = chr[0]
        loci = chr[1]
        if len(chrs) == 1:
            _gen_loci(stream, loci)
        else:
            for i in range(repeats):
                _gen_loci(stream, loci)

def generate_simcoal_from_template(model, chrs, params, out_dir = '.', tp_dir=None):
    '''
       Writes a complete SimCoal2 template file.

       This joins together get_demography_template and get_chr_template,
       which are feed into generate_model
       Please check the three functions for parameters (model from
         get_demography_template, chrs from get_chr_template and
         params from generate_model).
    '''
    stream = open(out_dir + sep + 'tmp.par', 'w')
    get_demography_template(stream, model, tp_dir)
    get_chr_template(stream, chrs)
    stream.close()
    #par_stream = open(out_dir + sep + 'tmp.par', 'r')
    #print par_stream.read()
    #par_stream.close()
    par_stream = open(out_dir + sep + 'tmp.par', 'r')
    generate_model(par_stream, model, params, out_dir = out_dir)
    par_stream.close()


