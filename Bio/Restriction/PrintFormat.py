#!/usr/bin/env python
#
#      Restriction Analysis Libraries.
#      Copyright (C) 2004. Frederic Sohm.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

import re
import itertools
from Bio.Restriction import RanaConfig as RanaConf

"""
Usage:

    PrintFormat allow to print the results from restriction analysis in 3
    different format.
    List, column or map.

    the easiest way to use it is:
    
    >>> from Bio.Restriction.PrintFormat import PrintFormat
    >>> from Bio.Restriction.Restriction import AllEnzymes
    >>> from Bio import Entrez
    >>> from Bio import SeqIO
    >>> handle = Entrez.efetch(db="nucleotide", rettype="fasta", id="SYNPBR322")
    >>> pBR322 = SeqIO.read(handle, "fasta")
    >>> handle.close()
    >>> dct = AllEnzymes.search(pBR322.seq)
    >>> new = PrintFormat()
    >>> new.print_that(dct, '\n my pBR322 analysis\n\n','\n no site :\n\n')

     my pBR322 analysis
     
    AasI       :  2169, 2582.
    AatII      :  4289.
    ...
            More enzymes.
    ...
    ZraI       :  4287.
    ZrmI       :  3847.
    
     no site:
     
    AarI      AatI      Acc65I    AcsI      AcvI      AdeI      AflII     AgeI    
    ...
            More enzymes.
    ...
    Vha464I   XapI      XbaI      XcmI      XhoI      XmaCI     XmaI      XmaJI     
    Zsp2I 

    >>> new.sequence = pBR322.seq
    >>> new.print_as("map")
    >>> new.print_that(dct)
    ...
    
    Some of the methods of PrintFormat are meant to be overriden by derived
    class.
"""

class PrintFormat(object):
    """PrintFormat allow the printing of results of restriction analysis."""

    ConsoleWidth = RanaConf.ConsoleWidth
    NameWidth    = RanaConf.NameWidth
    MaxSize      = RanaConf.MaxSize
    Cmodulo      = ConsoleWidth%NameWidth       
    PrefWidth    = ConsoleWidth - Cmodulo
    Indent       = RanaConf.Indent
    linesize     = PrefWidth - NameWidth

    def __init__(self):
        """PrintFormat() -> new PrintFormat Instance"""
        pass

    def print_as(self, what='list'):
        """PF.print_as([what='list']) -> print the results as specified.

        Valid format are:
            'list'      -> alphabetical order
            'number'    -> number of sites in the sequence
            'map'       -> a map representation of the sequence with the sites.

        If you want more flexibility over-ride the virtual method make_format.
        """
        if what == 'map':
            self.make_format = self._make_map
        elif what == 'number':
            self.make_format = self._make_number
        else:
            self.make_format = self._make_list
            
        return
            

    def print_that(self, dct, title='',  s1=''):
        """PF.print_that(dct, [title[, s1]]) -> Print dct nicely formatted.

        dct is a dictionary as returned by a RestrictionBatch.search()
        
        title is the title of the map.
        It must be a formated string, i.e. you must include the line break.
        
        s1 is the title separating the list of enzymes that have sites from
        those without sites.
        s1 must be a formatted string as well.

        The format of print_that is a list."""
        if not dct:
            dct = self.results
        ls, nc = [], []
        for k, v in dct.iteritems():
            if v:
                ls.append((k,v))
            else:
                nc.append(k)
        print self.make_format(ls, title, nc, s1)
        return
       
    def make_format(self, cut=[], title='', nc=[], s1=''):
        """PF.make_format(cut, nc, title, s) -> string

        Virtual method.
        Here to be pointed to one of the _make_* methods.
        You can as well create a new method and point make_format to it."""
        return self._make_list(cut,title, nc,s1)

###### _make_* methods to be used with the virtual method make_format

    def _make_list(self, ls,title, nc,s1):
        """PF._make_number(ls,title, nc,s1) -> string.

        return a string of form:
        
        title.

        enzyme1     :   position1, position2.
        enzyme2     :   position1, position2, position3.

        ls is a list of cutting enzymes.
        title is the title.
        nc is a list of non cutting enzymes.
        s1 is the sentence before the non cutting enzymes."""
        return self._make_list_only(ls, title) + self._make_nocut_only(nc, s1)

    def _make_map(self, ls,title, nc,s1):
        """PF._make_number(ls,title, nc,s1) -> string.

        return a string of form:
        
        title.

            enzyme1, position
            |
        AAAAAAAAAAAAAAAAAAAAA...
        |||||||||||||||||||||
        TTTTTTTTTTTTTTTTTTTTT...

        ls is a list of cutting enzymes.
        title is the title.
        nc is a list of non cutting enzymes.
        s1 is the sentence before the non cutting enzymes."""
        return self._make_map_only(ls, title) + self._make_nocut_only(nc, s1)

    def _make_number(self, ls,title, nc,s1):
        """PF._make_number(ls,title, nc,s1) -> string.

        title.
        
        enzyme which cut 1 time:
        
        enzyme1     :   position1.

        enzyme which cut 2 times:
        
        enzyme2     :   position1, position2.
        ...

        ls is a list of cutting enzymes.
        title is the title.
        nc is a list of non cutting enzymes.
        s1 is the sentence before the non cutting enzymes."""
        return self._make_number_only(ls, title)+self._make_nocut_only(nc,s1)
    
    def _make_nocut(self, ls,title, nc,s1):
        """PF._make_nocut(ls,title, nc,s1) -> string.

        return a formatted string of the non cutting enzymes.

        ls is a list of cutting enzymes -> will not be used.
        Here for compatibility with make_format.
        
        title is the title.
        nc is a list of non cutting enzymes.
        s1 is the sentence before the non cutting enzymes."""
        return title + self._make_nocut_only(nc, s1) 

    def _make_nocut_only(self, nc, s1, ls =[],title=''):
        """PF._make_nocut_only(nc, s1) -> string.

        return a formatted string of the non cutting enzymes.
        
        nc is a list of non cutting enzymes.
        s1 is the sentence before the non cutting enzymes.
        """
        if not nc:
            return s1
        nc.sort()
        st = ''
        stringsite = s1 or '\n   Enzymes which do not cut the sequence.\n\n'    
        Join = ''.join
        for key in nc:
            st = Join((st, str.ljust(str(key), self.NameWidth)))
            if len(st) > self.linesize:
                stringsite = Join((stringsite, st, '\n'))
                st = ''
        stringsite = Join((stringsite, st, '\n'))
        return stringsite
    
    def _make_list_only(self, ls, title, nc = [], s1 = ''):
        """PF._make_list_only(ls, title) -> string.

        return a string of form:
        
        title.

        enzyme1     :   position1, position2.
        enzyme2     :   position1, position2, position3.
        ...
        
        ls is a list of results.
        title is a string.
        Non cutting enzymes are not included."""
        if not ls:
            return title
        return self.__next_section(ls, title)

    def _make_number_only(self, ls, title, nc = [], s1 =''):
        """PF._make_number_only(ls, title) -> string.

        return a string of form:
        
        title.
        
        enzyme which cut 1 time:
        
        enzyme1     :   position1.

        enzyme which cut 2 times:
        
        enzyme2     :   position1, position2.
        ...
        
                
        ls is a list of results.
        title is a string.
        Non cutting enzymes are not included."""
        if not ls:
            return title
        ls.sort(lambda x,y : cmp(len(x[1]), len(y[1])))
        iterator = iter(ls)
        cur_len  = 1
        new_sect = []
        for name, sites in iterator:
            l = len(sites)
            if l > cur_len:
                title += "\n\nenzymes which cut %i times :\n\n"%cur_len
                title = self.__next_section(new_sect, title)
                new_sect, cur_len = [(name, sites)], l
                continue
            new_sect.append((name,sites))
        title += "\n\nenzymes which cut %i times :\n\n"%cur_len
        return self.__next_section(new_sect, title)
            
    def _make_map_only(self, ls, title, nc = [],  s1 = ''):
        """PF._make_map_only(ls, title) -> string.

        return a string of form:
        
        title.

            enzyme1, position
            |
        AAAAAAAAAAAAAAAAAAAAA...
        |||||||||||||||||||||
        TTTTTTTTTTTTTTTTTTTTT...
        
                
        ls is a list of results.
        title is a string.
        Non cutting enzymes are not included.
        """
        if not ls:
            return title
        resultKeys = [str(x) for x,y in ls]
        resultKeys.sort()
        map = title or ''
        enzymemap = {}
        for (enzyme, cut) in ls:
            for c in cut:
                if c in enzymemap:
                    enzymemap[c].append(str(enzyme))
                else:
                    enzymemap[c] = [str(enzyme)]
        mapping = enzymemap.keys()
        mapping.sort()
        cutloc = {}
        x, counter, length = 0, 0, len(self.sequence)
        for x in xrange(60, length, 60):
            counter = x - 60
            l=[]
            for key in mapping:
                if key <= x:
                    l.append(key)
                else:
                    cutloc[counter] = l
                    mapping = mapping[mapping.index(key):]
                    break
            cutloc[x] = l
        cutloc[x] = mapping
        sequence = self.sequence.tostring()
        revsequence = self.sequence.complement().tostring()
        a = '|'
        base, counter = 0, 0
        emptyline = ' ' * 60
        Join = ''.join
        for base in xrange(60, length, 60):
            counter = base - 60
            line = emptyline
            for key in cutloc[counter]:
                s = ''
                if key == base:
                    for n in enzymemap[key] : s = ' '.join((s,n))
                    l = line[0:59]
                    lineo = Join((l, str(key), s, '\n'))
                    line2 = Join((l, a, '\n'))
                    linetot = Join((lineo, line2))
                    map = Join((map, linetot))
                    break
                for n in enzymemap[key] : s = ' '.join((s,n))
                k = key%60
                lineo = Join((line[0:(k-1)], str(key), s, '\n'))
                line = Join((line[0:(k-1)], a, line[k:]))
                line2 = Join((line[0:(k-1)], a, line[k:], '\n'))
                linetot = Join((lineo,line2))
                map = Join((map,linetot))
            mapunit = '\n'.join((sequence[counter : base],a * 60,
                                 revsequence[counter : base],
                                 Join((str.ljust(str(counter+1), 15), ' '* 30,
                                        str.rjust(str(base), 15),'\n\n'))
                                 ))
            map = Join((map, mapunit)) 
        line = ' '* 60
        for key in cutloc[base]:
            s = ''
            if key == length:
                for n in enzymemap[key]:
                    s = Join((s,' ',n))
                l = line[0:(length-1)]
                lineo = Join((l,str(key),s,'\n'))
                line2 = Join((l,a,'\n'))
                linetot = Join((lineo, line2))
                map = Join((map, linetot))
                break
            for n in enzymemap[key] : s = Join((s,' ',n))
            k = key%60
            lineo = Join((line[0:(k-1)],str(key),s,'\n'))
            line = Join((line[0:(k-1)],a,line[k:]))
            line2 = Join((line[0:(k-1)],a,line[k:],'\n'))
            linetot = Join((lineo,line2))
            map = Join((map,linetot))
        mapunit = ''
        mapunit = Join((sequence[base : length], '\n'))
        mapunit = Join((mapunit, a * (length-base), '\n'))
        mapunit = Join((mapunit,revsequence[base:length], '\n'))
        mapunit = Join((mapunit, Join((str.ljust(str(base+1), 15), ' '*(
            length-base-30),str.rjust(str(length), 15),
                                       '\n\n'))))
        map = Join((map,mapunit))
        return map
    
###### private method to do lists:
    
    def __next_section(self, ls, into):
        """FP.__next_section(ls, into) -> string.

        ls is a list of tuple (string, [int, int]).
        into is a string to which the formatted ls will be added.

        Format ls as a string of lines:
        The form is:

        enzyme1     :   position1.
        enzyme2     :   position2, position3.

        then add the formatted ls to tot
        return tot."""
        ls.sort()
        indentation = '\n' + (self.NameWidth + self.Indent) * ' '
        linesize = self.linesize - self.MaxSize
        pat = re.compile("([\w,\s()]){1,%i}[,\.]"%linesize)
        several, Join = '', ''.join
        for name, sites in ls:
            stringsite = ''
            l = Join((', '.join([str(site) for site in sites]), '.'))
            if len(l) > linesize:
                #
                #   cut where appropriate and add the indentation
                #
                l = [x.group() for x in re.finditer(pat, l)]
                stringsite = indentation.join(l) 
            else:
                stringsite = l    
            into = Join((into,
                         str(name).ljust(self.NameWidth),' :  ',stringsite,'\n'))
        return into
