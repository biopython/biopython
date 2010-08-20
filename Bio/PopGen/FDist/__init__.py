# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""
This module provides code to work with FDist.

See http://www.rubic.rdg.ac.uk/~mab/software.html .

Classes:
Record           Holds FDist data.

Functions:
read             Parses a FDist record (file) into a Record object.


Deprecated classes:
RecordParser     Parses a FDist record (file) into a Record object.

_Scanner         Scans a FDist record.
_RecordConsumer  Consumes FDist data to a Record object.


"""




def read(handle):
    """Parses FDist data into a Record object.

       handle is a file-like object that contains a FDist record.
    """
    record = Record()
    record.data_org = int(str(handle.next()).rstrip())
    record.num_pops = int(str(handle.next()).rstrip())
    record.num_loci = int(str(handle.next()).rstrip())
    for i in range(record.num_loci):
        handle.next()
        num_alleles = int(str(handle.next()).rstrip())
        pops_data = []
        if record.data_org==0:
            for j in range(record.num_pops):
                line_comp = str(handle.next()).rstrip().split(' ')
                pop_dist = map(lambda x: int(x), line_comp)
                pops_data.append(pop_dist)
        else:
            raise NotImplementedError('1/alleles by rows not implemented')
        record.loci_data.append((num_alleles, pops_data))
    return record


class Record:
    """Holds information from a FDist record.

    Members:
    data_org    Data organization (0 pops by rows, 1 alleles by rows).
                The Record will behave as if data was 0 (converting if needed)
    
    num_pops       Number of populations
    
    num_loci       Number of loci
    
    loci_data      Loci data
    
    loci_data is a list, where each element represents a locus. Each element
    is a tuple, the first element is the number of alleles, the second
    element a list. Each element of the list is the count of each allele
    per population.
    
    """
    def __init__(self):
        self.data_org    = 0
        self.num_pops    = 0
        self.num_loci    = 0
        self.loci_data   = []
        
    def __str__(self):
        rep  = ['0\n'] #We only export in 0 format, even if originally was 1
        rep.append(str(self.num_pops) + '\n')
        rep.append(str(self.num_loci) + '\n')
        rep.append('\n')
        for locus_data in self.loci_data:
            num_alleles, pops_data = locus_data
            rep.append(str(num_alleles) + '\n')
            for pop_data in pops_data:
                for allele_count in pop_data:
                    rep.append(str(allele_count) + ' ')
                rep.append('\n')
            rep.append('\n')
        return "".join(rep)
    
# Everything below is deprecated

from Bio.ParserSupport import *

class RecordParser(AbstractParser):
    """Parses FDist data into a Record object (DEPRECATED).

    This class is DEPRECATED; please use the read() function in this
    module instead.
    """
    def __init__(self):
        import warnings
        from Bio import BiopythonDeprecationWarning
        warnings.warn("Bio.PopGen.FDist.RecordParser is deprecated. Please use the read() function in this module instead.",
              BiopythonDeprecationWarning)
        self._scanner = _Scanner()
        self._consumer = _RecordConsumer()

    def parse(self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class _Scanner:
    """Scans a FDist record (DEPRECATED).
    
    There is only one record per file.

    This class is DEPRECATED; please use the read() function in this
    module instead.
    """

    def __init__(self):
        import warnings
        from Bio import BiopythonDeprecationWarning
        warnings.warn("Bio.PopGen.FDist._Scanner is deprecated. Please use the read() function in this module instead.",
              BiopythonDeprecationWarning)

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in a FDist unit record for scanning.  handle is a file-like
        object that contains a FDist record.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """

        consumer.start_record()
        self.num_pops    = None
        self.num_loci    = None
        self.loci_data   = []
        
        data_org = int(handle.readline().rstrip())
        consumer.data_org(data_org)
        num_pops = int(handle.readline().rstrip())
        consumer.num_pops(num_pops)
        num_loci = int(handle.readline().rstrip())
        consumer.num_loci(num_loci)
        for i in range(num_loci):
            handle.readline()
            num_alleles = int(handle.readline().rstrip())
            pops_data = []
            if data_org==0:
                for j in range(num_pops):
                    line_comp = handle.readline().rstrip().split(' ')
                    pop_dist = map(lambda x: int(x), line_comp)
                    pops_data.append(pop_dist)
            else:
                raise NotImplementedError('1/alleles by rows not implemented')
            consumer.new_locus(num_alleles, pops_data)
        consumer.end_record()

class _RecordConsumer(AbstractConsumer):
    """Consumer that converts a FDist record to a Record object (DEPRECATED).

    Members:
    data    Record with FDist data.

    This class is DEPRECATED; please use the read() function in this module
    instead.
    """
    def __init__(self):
        import warnings
        from Bio import BiopythonDeprecationWarning
        warnings.warn("Bio.PopGen.FDist._RecordConsumer is deprecated. Please use the read() function in this module instead.",
              BiopythonDeprecationWarning)
        self.data = None

    def start_record(self):
        self.data = Record()

    def end_record(self):
        pass

    def data_org(self, data_org):
        self.data.data_org = data_org
        
    def num_pops(self, num_pops):
        self.data.num_pops = num_pops

    def num_loci(self, num_loci):
        self.data.num_loci = num_loci

    def new_locus(self, num_alleles, pop_data):
        self.data.loci_data.append((num_alleles, pop_data))
