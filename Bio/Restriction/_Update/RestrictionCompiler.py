#!/usr/bin/env python
#
#      Restriction Analysis Libraries.
#      Copyright (C) 2004. Frederic Sohm.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
#   this script is used to produce the dictionary which will contains the data
#   about the restriction enzymes from the Emboss/Rebase data files
#   namely
#   emboss_e.### (description of the sites),
#   emboss_r.### (origin, methylation, references)
#   emboss_s.### (suppliers)
#   where ### is a number of three digits : 1 for the year two for the month
#
#   very dirty implementation but it does the job, so...
#   Not very quick either but you are not supposed to use it frequently.
#
#   The results are stored in
#   path/to/site-packages/Bio/Restriction/Restriction_Dictionary.py
#   the file contains two dictionary:
#   'rest_dict' which contains the data for the enzymes
#   and
#   'suppliers' which map the name of the suppliers to their abbreviation.
#

"""Convert a serie of Rebase files into a Restriction_Dictionary.py module.

The Rebase files are in the emboss format:

    emboss_e.###    -> contains informations about the restriction sites.
    emboss_r.###    -> contains general informations about the enzymes.
    emboss_s.###    -> contains informations about the suppliers.
    
### is a 3 digit number. The first digit is the year and the two last the month.
"""

import sre
import os
import itertools
import time
import sys
import site
import shutil

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import Bio.Restriction.Restriction
from Bio.Restriction.Restriction import AbstractCut, RestrictionType, NoCut, OneCut,\
TwoCuts, Meth_Dep, Meth_Undep, Palindromic, NonPalindromic, Unknown, Blunt,\
Ov5, Ov3, NotDefined, Defined, Ambiguous, Commercially_available, Not_available 

import Bio.Restriction.RanaConfig as config
from Bio.Restriction._Update.Update import RebaseUpdate
from Bio.Restriction.Restriction import *
from Bio.Restriction.DNAUtils import antiparallel

DNA=Seq
dna_alphabet = {'A':'A', 'C':'C', 'G':'G', 'T':'T',
                'R':'AG', 'Y':'CT', 'W':'AT', 'S':'CG', 'M':'AC', 'K':'GT',
                'H':'ACT', 'B':'CGT', 'V':'ACG', 'D':'AGT',
                'N':'ACGT',
                'a': 'a', 'c': 'c', 'g': 'g', 't': 't',
                'r':'ag', 'y':'ct', 'w':'at', 's':'cg', 'm':'ac', 'k':'gt',
                'h':'act', 'b':'cgt', 'v':'acg', 'd':'agt',
                'n':'acgt'}


complement_alphabet = {'A':'T', 'T':'A', 'C':'G', 'G':'C','R':'Y', 'Y':'R',
                       'W':'W', 'S':'S', 'M':'K', 'K':'M', 'H':'D', 'D':'H',
                       'B':'V', 'V':'B', 'N':'N','a':'t', 'c':'g', 'g':'c',
                       't':'a', 'r':'y', 'y':'r', 'w':'w', 's':'s','m':'k',
                       'k':'m', 'h':'d', 'd':'h', 'b':'v', 'v':'b', 'n':'n'}
enzymedict = {}
suppliersdict = {}
classdict = {}
typedict = {}


class OverhangError(ValueError):
    """Exception for dealing with overhang."""
    pass
          
def BaseExpand(base):
    """BaseExpand(base) -> string.

    given a degenerated base, returns its meaning in IUPAC alphabet.

    i.e:
        b= 'A' -> 'A'
        b= 'N' -> 'ACGT'
        etc..."""
    base = base.upper()
    return dna_alphabet[base]

def regex(site):
    """regex(site) -> string.

    Construct a regular expression from a DNA sequence.
    i.e.:
        site = 'ABCGN'   -> 'A[CGT]CG.'"""
    reg_ex = site
    for base in reg_ex:
        if base in ('A', 'T', 'C', 'G', 'a', 'c', 'g', 't'):
            pass
        if base in ('N', 'n'):
            reg_ex = '.'.join(reg_ex.split('N'))
            reg_ex = '.'.join(reg_ex.split('n'))
        if base in ('R', 'Y', 'W', 'M', 'S', 'K', 'H', 'D', 'B', 'V'):
            expand = '['+ str(BaseExpand(base))+']'
            reg_ex = expand.join(reg_ex.split(base))
    return reg_ex

def Antiparallel(sequence):
    """Antiparallel(sequence) -> string.

    returns a string which represents the reverse complementary strand of
    a DNA sequence."""
    return antiparallel(sequence.tostring())

def is_palindrom(sequence):
    """is_palindrom(sequence) -> bool.

    True is the sequence is a palindrom.
    sequence is a DNA object."""
    return sequence == DNA(Antiparallel(sequence))

def LocalTime():
    """LocalTime() -> string.

    LocalTime calculate the extension for emboss file for the current year and
    month."""
    t = time.gmtime()
    year = str(t.tm_year)[-1]
    month = str(t.tm_mon)
    if len(month) == 1 : month = '0'+month
    return year+month
                

class newenzyme(object):
    """construct the attributes of the enzyme corresponding to 'name'."""
    def __init__(cls, name):
        cls.opt_temp = 37
        cls.inact_temp = 65
        cls.substrat = 'DNA'
        target = enzymedict[name]
        cls.site = target[0]
        cls.size = target[1]
        cls.suppl = tuple(target[9])
        cls.freq = target[11]
        cls.ovhg = target[13]
        cls.ovhgseq = target[14]
        cls.bases = ()
        #
        #   Is the site palindromic?
        #   Important for the way the DNA is search for the site.
        #   Palindromic sites needs to be looked for only over 1 strand.
        #   Non Palindromic needs to be search for on the reverse complement
        #   as well.
        #
        if target[10] : cls.bases += ('Palindromic',)
        else : cls.bases += ('NonPalindromic',)
        #
        #   Number of cut the enzyme produce.
        #   0 => unknown, the enzyme has not been fully characterised.
        #   2 => 1 cut, (because one cut is realised by cutting 2 strands
        #   4 => 2 cuts, same logic.
        #   A little bit confusing but it is the way EMBOSS/Rebase works.
        #
        if not target[2]:
            #
            #   => undefined enzymes, nothing to be done.
            #
            cls.bases += ('NoCut','Unknown', 'NotDefined')
            cls.fst5 = None
            cls.fst3 = None
            cls.scd5 = None
            cls.scd3 = None
            cls.ovhg = None
            cls.ovhgseq = None
        else:
            #
            #   we will need to calculate the overhang.
            #
            if target[2] == 2:
                cls.bases += ('OneCut',)
                cls.fst5 = target[4]
                cls.fst3 = target[5]
                cls.scd5 = None
                cls.scd3 = None
            else:
                cls.bases += ('TwoCuts',)
                cls.fst5 = target[4]
                cls.fst3 = target[5]
                cls.scd5 = target[6]
                cls.scd3 = target[7]
            #
            #   Now, prepare the overhangs which will be added to the DNA 
            #   after the cut.
            #   Undefined enzymes will not be allowed to catalyse,
            #   they are not available commercially anyway.
            #   I assumed that if an enzyme cut twice the overhang will be of
            #   the same kind. The only exception is HaeIV. I do not deal
            #   with that at the moment (ie I don't include it,
            #   need to be fixed).
            #   They generally cut outside their recognition site and
            #   therefore the overhang is undetermined and dependent of
            #   the DNA sequence upon which the enzyme act.
            #
            if target[3]:
                #
                #   rebase field for blunt: blunt == 1, other == 0.
                #   The enzyme is blunt. No overhang.
                #
                cls.bases += ('Blunt', 'Defined')
                cls.ovhg = 0
            elif isinstance(cls.ovhg, int):
                #
                #   => overhang is sequence dependent
                #
                if cls.ovhg > 0:
                    #
                    #   3' overhang, ambiguous site (outside recognition site
                    #   or site containing ambiguous bases (N, W, R,...)
                    #
                    cls.bases += ('Ov3', 'Ambiguous')
                elif cls.ovhg < 0:
                    #
                    #   5' overhang, ambiguous site (outside recognition site
                    #   or site containing ambiguous bases (N, W, R,...)
                    #
                    cls.bases += ('Ov5', 'Ambiguous')
            else:
                #
                #   cls.ovhg is a string => overhang is constant
                #
                if cls.fst5 - (cls.fst3 + cls.size) < 0:
                    cls.bases += ('Ov5', 'Defined')
                    cls.ovhg = - len(cls.ovhg)
                else:
                    cls.bases += ('Ov3', 'Defined')
                    cls.ovhg = + len(cls.ovhg)
        #
        #   Next class : sensibility to methylation.
        #   Set by EmbossMixer from emboss_r.txt file
        #   Not really methylation dependent at the moment, stands rather for
        #   'is the site methylable?'.
        #   Proper methylation sensibility has yet to be implemented.
        #   But the class is there for further development.
        #
        if target[8]:
            cls.bases += ('Meth_Dep', )
            cls.compsite = target[12]
        else:
            cls.bases += ('Meth_Undep',)
            cls.compsite = target[12]
        #
        #   Next class will allow to select enzymes in function of their
        #   suppliers. Not essential but can be useful.
        #
        if cls.suppl:
            cls.bases += ('Commercially_available', )
        else:
            cls.bases += ('Not_available', )
        cls.bases += ('AbstractCut', 'RestrictionType')
        cls.__name__ = name
        cls.results = None
        cls.dna = None
        cls.__bases__ = cls.bases
        cls.charac = (cls.fst5, cls.fst3, cls.scd5, cls.scd3, cls.site)
        if not target[2] and cls.suppl:
            supp = ', '.join([suppliersdict[s][0] for s in cls.suppl])
            print 'WARNING : It seems that %s is both commercially available\
            \n\tand its characteristics are unknown. \
            \n\tThis seems counter-intuitive.\
            \n\tThere is certainly an error either in ranacompiler or\
            \n\tin this REBASE release.\
            \n\tThe supplier is : %s.' % (name, supp)
        return


class TypeCompiler(object):
    
    """Build the different types possible for Restriction Enzymes"""

    def __init__(self):
        """TypeCompiler() -> new TypeCompiler instance."""
        pass

    def buildtype(self):
        """TC.buildtype() -> generator.

        build the new types that will be needed for constructing the
        restriction enzymes."""
        baT = (AbstractCut, RestrictionType)
        cuT = (NoCut, OneCut, TwoCuts)
        meT = (Meth_Dep, Meth_Undep)
        paT = (Palindromic, NonPalindromic)
        ovT = (Unknown, Blunt, Ov5, Ov3)
        deT = (NotDefined, Defined, Ambiguous)
        coT = (Commercially_available, Not_available)
        All = (baT, cuT, meT, paT, ovT, deT, coT)
        #
        #   Now build the types. Only the most obvious are left out.
        #   Modified even the most obvious are not so obvious.
        #   emboss_*.403 AspCNI is unknown and commercially available.
        #   So now do not remove the most obvious.
        #
        types = [(p,c,o,d,m,co,baT[0],baT[1])
                 for p in paT for c in cuT for o in ovT
                 for d in deT for m in meT for co in coT]
        n= 1
        for ty in types:
            dct = {}
            for t in ty:
                dct.update(t.__dict__)
                #
                #   here we need to customize the dictionary.
                #   i.e. types deriving from OneCut have always scd5 and scd3
                #   equal to None. No need therefore to store that in a specific
                #   enzyme of this type. but it then need to be in the type.
                #
                dct['results'] = []
                dct['substrat'] = 'DNA'
                dct['dna'] = None
                if t == NoCut:
                    dct.update({'fst5':None,'fst3':None,
                                'scd5':None,'scd3':None,
                                'ovhg':None,'ovhgseq':None})
                elif t == OneCut:
                    dct.update({'scd5':None, 'scd3':None})
            class klass(type):
                def __new__(cls):
                    return type.__new__(cls, 'type%i'%n,ty,dct)
                def __init__(cls):
                    super(klass, cls).__init__('type%i'%n,ty,dct)
            yield klass()
            n+=1

start = '\n\
#!/usr/bin/env python\n\
#\n\
#      Restriction Analysis Libraries.\n\
#      Copyright (C) 2004. Frederic Sohm.\n\
#\n\
# This code is part of the Biopython distribution and governed by its\n\
# license.  Please see the LICENSE file that should have been included\n\
# as part of this package.\n\
#\n\
# This file is automatically generated - do not edit it by hand! Instead,\n\
# use the tool Scripts/Restriction/ranacompiler.py which in turn uses\n\
# Bio/Restriction/_Update/RestrictionCompiler.py\n\
#\n\
# The following dictionaries used to be defined in one go, but that does\n\
# not work on Jython due to JVM limitations. Therefore we break this up\n\
# into steps, using temporary functions to avoid the JVM limits.\n\
\n\n'

class DictionaryBuilder(object):

    def __init__(self, e_mail='', ftp_proxy=''):
        """DictionaryBuilder([e_mail[, ftp_proxy]) -> DictionaryBuilder instance.

        If the emboss files used for the construction need to be updated this
        class will download them if the ftp connection is correctly set.
        either in RanaConfig.py or given at run time.
        
        e_mail is the e-mail address used as password for the anonymous
        ftp connection.

        proxy is the ftp_proxy to use if any."""
        self.rebase_pass = e_mail or config.Rebase_password
        self.proxy = ftp_proxy or config.ftp_proxy
    
    def build_dict(self):
        """DB.build_dict() -> None.

        Construct the dictionary and build the files containing the new
        dictionaries."""
        #
        #   first parse the emboss files.
        #
        emboss_e, emboss_r, emboss_s = self.lastrebasefile()
        #
        #   the results will be stored into enzymedict.
        #
        self.information_mixer(emboss_r, emboss_e, emboss_s)
        emboss_r.close()
        emboss_e.close()
        emboss_s.close()
        #
        #   we build all the possible type 
        #
        tdct = {}
        for klass in TypeCompiler().buildtype():
            exec klass.__name__ +'= klass'
            exec "tdct['"+klass.__name__+"'] = klass"

        #
        #   Now we build the enzymes from enzymedict
        #   and store them in a dictionary.
        #   The type we will need will also be stored.
        #

        for name in enzymedict:
            #
            #   the class attributes first:
            #
            cls = newenzyme(name)
            #
            #   Now select the right type for the enzyme.
            #
            bases = cls.bases
            clsbases = tuple([eval(x) for x in bases])
            typestuff = ''
            for n, t in tdct.iteritems():
                #
                #   if the bases are the same. it is the right type.
                #   create the enzyme and remember the type
                #
                if t.__bases__ == clsbases:
                    typestuff = t
                    typename = t.__name__
                continue
            #
            #   now we build the dictionaries.
            #
            dct = dict(cls.__dict__)
            del dct['bases']
            del dct['__bases__']
            del dct['__name__']# no need to keep that, it's already in the type.
            classdict[name] = dct
           
            commonattr = ['fst5', 'fst3', 'scd5', 'scd3', 'substrat',
                          'ovhg', 'ovhgseq','results', 'dna']
            if typename in typedict:
                typedict[typename][1].append(name)
            else:
                enzlst= []
                tydct = dict(typestuff.__dict__)
                tydct = dict([(k,v) for k,v in tydct.iteritems() if k in commonattr])
                enzlst.append(name)
                typedict[typename] = (bases, enzlst)
            for letter in cls.__dict__['suppl']:
                supplier = suppliersdict[letter]
                suppliersdict[letter][1].append(name)
        if not classdict or not suppliersdict or not typedict:
            print 'One of the new dictionaries is empty.'
            print 'Check the integrity of the emboss file before continuing.'
            print 'Update aborted.'
            sys.exit()
        #
        #   How many enzymes this time?
        #
        print '\nThe new database contains %i enzymes.\n' % len(classdict)
        #
        #   the dictionaries are done. Build the file 
        #
        #update = config.updatefolder
        
        update = os.getcwd()
        results = open(os.path.join(update, 'Restriction_Dictionary.py'), 'w')
        print 'Writing the dictionary containing the new Restriction classes.\t',
        results.write(start)
        results.write('rest_dict = {}\n')
        for name in sorted(classdict) :
            results.write("def _temp():\n")
            results.write("    return {\n")
            for key, value in classdict[name].iteritems() :
                results.write("        %s : %s,\n" % (repr(key), repr(value)))
            results.write("    }\n")
            results.write("rest_dict[%s] = _temp()\n" % repr(name))
            results.write("\n")
        print 'OK.\n'
        print 'Writing the dictionary containing the suppliers datas.\t\t',
        results.write('suppliers = {}\n')
        for name in sorted(suppliersdict) :
            results.write("def _temp():\n")
            results.write("    return (\n")
            for value in suppliersdict[name] :
                results.write("        %s,\n" % repr(value))
            results.write("    )\n")
            results.write("suppliers[%s] = _temp()\n" % repr(name))
            results.write("\n")
        print 'OK.\n'
        print 'Writing the dictionary containing the Restriction types.\t',
        results.write('typedict = {}\n')
        for name in sorted(typedict) :
            results.write("def _temp():\n")
            results.write("    return (\n")
            for value in typedict[name] :
                results.write("        %s,\n" % repr(value))
            results.write("    )\n")
            results.write("typedict[%s] = _temp()\n" % repr(name))
            results.write("\n")
        #I had wanted to do "del _temp" at each stage (just for clarity), but
        #that pushed the code size just over the Jython JVM limit. We include
        #one the final "del _temp" to clean up the namespace.
        results.write("del _temp\n")
        results.write("\n")
        print 'OK.\n'
        results.close()
        return

    def install_dict(self):
        """DB.install_dict() -> None.

        Install the newly created dictionary in the site-packages folder.

        May need super user privilege on some architectures."""
        print '\n ' +'*'*78 + ' \n'
        print '\n\t\tInstalling Restriction_Dictionary.py'
        try:
            import Bio.Restriction.Restriction_Dictionary as rd
        except ImportError:
            print '\
            \n Unable to locate the previous Restriction_Dictionary.py module\
            \n Aborting installation.'
            sys.exit()
        #
        #   first save the old file in Updates
        #
        old = os.path.join(os.path.split(rd.__file__)[0],
                           'Restriction_Dictionary.py')
        #update_folder = config.updatefolder
        update_folder = os.getcwd()
        shutil.copyfile(old, os.path.join(update_folder,
                                          'Restriction_Dictionary.old'))
        #
        #   Now test and install.
        #
        new = os.path.join(update_folder, 'Restriction_Dictionary.py')
        try:
            execfile(new)
            print '\
            \n\tThe new file seems ok. Proceeding with the installation.'   
        except SyntaxError:
            print '\
            \n The new dictionary file is corrupted. Aborting the installation.'
            return
        try:
            shutil.copyfile(new, old)
            print'\n\t Everything ok. If you need it a version of the old\
            \n\t dictionary have been saved in the Updates folder under\
            \n\t the name Restriction_Dictionary.old.'
            print '\n ' +'*'*78 + ' \n'
        except IOError:
            print '\n ' +'*'*78 + ' \n'
            print '\
            \n\t WARNING : Impossible to install the new dictionary.\
            \n\t Are you sure you have write permission to the folder :\n\
            \n\t %s ?\n\n' % os.path.split(old)[0]
            return self.no_install()
        return

    def no_install(self):
        """BD.no_install() -> None.

        build the new dictionary but do not install the dictionary."""
        print '\n ' +'*'*78 + '\n'
        #update = config.updatefolder
        try:
            import Bio.Restriction.Restriction_Dictionary as rd
        except ImportError:
            print '\
            \n Unable to locate the previous Restriction_Dictionary.py module\
            \n Aborting installation.'
            sys.exit()
        #
        #   first save the old file in Updates
        #
        old = os.path.join(os.path.split(rd.__file__)[0],
                           'Restriction_Dictionary.py')
        update = os.getcwd()
        shutil.copyfile(old, os.path.join(update, 'Restriction_Dictionary.old'))
        places = update, os.path.split(Bio.Restriction.Restriction.__file__)[0]
        print "\t\tCompilation of the new dictionary : OK.\
        \n\t\tInstallation : No.\n\
        \n You will find the newly created 'Restriction_Dictionary.py' file\
        \n in the folder : \n\
        \n\t%s\n\
        \n Make a copy of 'Restriction_Dictionary.py' and place it with \
        \n the other Restriction libraries.\n\
        \n note : \
        \n This folder should be :\n\
        \n\t%s\n" % places
        print '\n ' +'*'*78 + '\n'
        return
        

    def lastrebasefile(self):
        """BD.lastrebasefile() -> None.

        Check the emboss files are up to date and download them if they are not.
        """
        embossnames = ('emboss_e', 'emboss_r', 'emboss_s')
        #
        #   first check if we have the last update:
        #
        emboss_now = ['.'.join((x,LocalTime())) for x in embossnames]
        update_needed = False
        #dircontent = os.listdir(config.Rebase) #    local database content
        dircontent = os.listdir(os.getcwd())
        base = os.getcwd() # added for biopython current directory
        for name in emboss_now:
            if name in dircontent:
                pass
            else:
                update_needed = True
                
        if not update_needed:
            #
            #   nothing to be done
            #
            print '\n Using the files : %s'% ', '.join(emboss_now)
            return tuple([open(os.path.join(base, n)) for n in emboss_now])
        else:
            #
            #   may be download the files.
            #
            print '\n The rebase files are more than one month old.\
            \n Would you like to update them before proceeding?(y/n)'
            r = raw_input(' update [n] >>> ')
            if r in ['y', 'yes', 'Y', 'Yes']:
                updt = RebaseUpdate(self.rebase_pass, self.proxy)
                updt.openRebase()
                updt.getfiles()
                updt.close()
                print '\n Update complete. Creating the dictionaries.\n'
                print '\n Using the files : %s'% ', '.join(emboss_now)
                return tuple([open(os.path.join(base, n)) for n in emboss_now])
            else:
                #
                #   we will use the last files found without updating.
                #   But first we check we have some file to use.
                #
                class NotFoundError(Exception):
                    pass
                for name in embossnames:
                    try:
                        for file in dircontent:
                            if file.startswith(name):
                                break
                        else:
                            pass
                        raise NotFoundError
                    except NotFoundError:
                        print "\nNo %s file found. Upgrade is impossible.\n"%name
                        sys.exit()
                    continue
                pass
        #
        #   now find the last file.
        #
        last = [0]
        for file in dircontent:
            fs = file.split('.')
            try:
                if fs[0] in embossnames and int(fs[1]) > int(last[-1]):
                    if last[0] : last.append(fs[1])
                    else : last[0] = fs[1]
                else:
                    continue
            except ValueError:
                continue
        last.sort()
        last = last[::-1]
        if int(last[-1]) < 100 : last[0], last[-1] = last[-1], last[0]
        for number in last:
            files = [(name, name+'.%s'%number) for name in embossnames]
            strmess = '\nLast EMBOSS files found are :\n'
            try:
                for name,file in files:
                    if os.path.isfile(os.path.join(base, file)):
                        strmess += '\t%s.\n'%file
                    else:
                        raise ValueError
                print strmess
                emboss_e = open(os.path.join(base, 'emboss_e.%s'%number),'r')
                emboss_r = open(os.path.join(base, 'emboss_r.%s'%number),'r')
                emboss_s = open(os.path.join(base, 'emboss_s.%s'%number),'r')
                return emboss_e, emboss_r, emboss_s
            except ValueError:
                continue

    def parseline(self, line):
        line = [line[0]]+[line[1].upper()]+[int(i) for i in line[2:9]]+line[9:]
        name = line[0].replace("-","_")
        site = line[1]          #   sequence of the recognition site
        dna = DNA(site)  
        size = line[2]          #   size of the recognition site
        #
        #   Calculate the overhang.
        #
        fst5 = line[5]  #   first site sense strand 
        fst3 = line[6]  #   first site antisense strand
        scd5 = line[7]  #   second site sense strand
        scd3 = line[8]  #   second site antisense strand
        
        #
        #   the overhang is the difference between the two cut
        #
        ovhg1 = fst5 - fst3
        ovhg2 = scd5 - scd3
        
        #
        #   0 has the meaning 'do not cut' in rebase. So we get short of 1
        #   for the negative numbers so we add 1 to negative sites for now.
        #   We will deal with the record later.
        #
        
        if fst5 < 0 : fst5 += 1
        if fst3 < 0 : fst3 += 1
        if scd5 < 0 : scd5 += 1
        if scd3 < 0 : scd3 += 1

        if ovhg2 != 0 and ovhg1 != ovhg2:
            #
            #   different length of the overhang of the first and second cut
            #   it's a pain to deal with and at the moment it concerns only
            #   one enzyme which is not commercially available (HaeIV).
            #   So we don't deal with it but we check the progression
            #   of the affair.
            #   Should HaeIV become commercially available or other similar
            #   new enzymes be added, this might be modified.
            #
            print '\
            \nWARNING : %s cut twice with different overhang length each time.\
            \n\tUnable to deal with this behaviour. \
            \n\tThis enzyme will not be included in the database. Sorry.' %name
            print '\tChecking :',
            raise OverhangError
        if 0 <= fst5 <= size and 0 <= fst3 <= size:
            #
            # cut inside recognition site
            #
            if fst5 < fst3:
                #
                #  5' overhang
                #
                ovhg1 = ovhgseq = site[fst5:fst3]       
            elif fst5 > fst3:
                #
                #  3' overhang
                #
                ovhg1 = ovhgseq = site[fst3:fst5]  
            else:
                #
                #  blunt
                #
                ovhg1 = ovhgseq = ''            
            for base in 'NRYWMSKHDBV':
                if base in ovhg1:
                    #
                    #   site and overhang degenerated
                    #
                    ovhgseq = ovhg1
                    if fst5 < fst3 :  ovhg1 = - len(ovhg1)
                    else : ovhg1 = len(ovhg1)
                    break
                else:
                    continue
        elif 0 <= fst5 <= size:
            #
            #   5' cut inside the site 3' outside
            #
            if fst5 < fst3:
                #
                #   3' cut after the site
                #
                ovhgseq = site[fst5:] + (fst3 - size) * 'N' 
            elif fst5 > fst3:
                #
                #   3' cut before the site
                #
                ovhgseq = abs(fst3) * 'N' + site[:fst5] 
            else:
                #
                #   blunt outside
                #
                ovhg1 = ovhgseq = '' 
        elif 0 <= fst3 <= size:
            #
            #   3' cut inside the site, 5' outside
            #
            if fst5 < fst3:
                #
                #   5' cut before the site
                #
                ovhgseq = abs(fst5) * 'N' + site[:fst3]
            elif fst5 > fst3:
                #
                #   5' cut after the site
                #
                ovhgseq = site[fst3:] + (fst5 - size) * 'N'
            else:
                #
                #   should not happend
                #
                raise ValueError('Error in #1')
        elif fst3 < 0 and size < fst5:
            #
            #   3' overhang. site is included.
            #
            ovhgseq = abs(fst3)*'N' + site + (fst5-size)*'N'
        elif fst5 < 0 and size <fst3:
            #
            #   5' overhang. site is included.
            #
            ovhgseq = abs(fst5)*'N' + site + (fst3-size)*'N'
        else:
            #
            #   5' and  3' outside of the site
            #
            ovhgseq = 'N' * abs(ovhg1)
        #
        #   Now line[5] to [8] are the location of the cut but we have to
        #   deal with the weird mathematics of biologists.
        #
        #   EMBOSS sequence numbering give:
        #                 DNA = 'a c g t A C G T'
        #                             -1 1 2 3 4
        #
        #   Biologists do not know about 0. Too much use of latin certainly.
        #
        #   To compensate, we add 1 to the positions if they are negative.
        #   No need to modify 0 as it means no cut and will not been used.
        #   Positive numbers should be ok since our sequence starts 1.
        #
        #   Moreover line[6] and line[8] represent cut on the reverse strand.
        #   They will be used for non palindromic sites and sre.finditer
        #   will detect the site in inverse orientation so we need to add the
        #   length of the site to compensate (+1 if they are negative).
        #
        for x in (5, 7):
            if line[x] < 0 : line[x] += 1
        for x in (6, 8):
            if line[x] > 0 : line[x] -= size
            elif line[x] < 0 : line[x] = line[x] - size + 1
        #
        #   now is the site palindromic?
        #   produce the regular expression which correspond to the site.
        #   tag of the regex will be the name of the enzyme for palindromic
        #   enzymesband two tags for the other, the name for the sense sequence
        #   and the name with '_as' at the end for the antisense sequence.
        #
        rg = ''
        if is_palindrom(dna):
            line.append(True)
            rg = ''.join(['(?P<', name, '>', regex(site.upper()), ')'])
        else:
            line.append(False)
            sense = ''.join(['(?P<', name, '>', regex(site.upper()), ')'])
            antisense = ''.join(['(?P<', name, '_as>',
                                 regex(Antiparallel(dna)), ')'])
            rg = sense + '|' + antisense
        #
        #   exact frequency of the site. (ie freq(N) == 1, ...)
        #
        f = [4/len(dna_alphabet[l]) for l in site.upper()]
        freq = reduce(lambda x, y : x*y, f)
        line.append(freq)
        #
        #   append regex and ovhg1, they have not been appended before not to
        #   break the factory class. simply to leazy to make the changes there.
        #
        line.append(rg)
        line.append(ovhg1)
        line.append(ovhgseq)
        return line

    def removestart(self, file):
        #
        #   remove the heading of the file.
        #
        return [l for l in itertools.dropwhile(lambda l:l.startswith('#'),file)]

    def getblock(self, file, index):
        #
        #   emboss_r.txt, separation between blocks is //
        #
        take = itertools.takewhile
        block = [l for l in take(lambda l :not l.startswith('//'),file[index:])]
        index += len(block)+1
        return block, index

    def get(self, block):
        #
        #   take what we want from the block.
        #   Each block correspond to one enzyme.
        #   block[0] => enzyme name
        #   block[3] => methylation (position and type)
        #   block[5] => suppliers (as a string of single letter)
        #
        bl3 = block[3].strip() 
        if not bl3 : bl3 = False #  site is not methylable
        return (block[0].strip(), bl3, block[5].strip())

    def information_mixer(self, file1, file2, file3):
        #
        #   Mix all the information from the 3 files and produce a coherent
        #   restriction record.
        #     
        methfile = self.removestart(file1)
        sitefile = self.removestart(file2)
        supplier = self.removestart(file3)
        
        i1, i2= 0, 0
        try:
            while True:
                block, i1 = self.getblock(methfile, i1)
                bl = self.get(block)
                line = (sitefile[i2].strip()).split()
                name = line[0]
                if name == bl[0]:
                    line.append(bl[1])  #   -> methylation
                    line.append(bl[2])  #   -> suppliers
                else:
                    bl = self.get(oldblock)
                    if line[0] == bl[0]:
                        line.append(bl[1])
                        line.append(bl[2])
                        i2 += 1
                    else:
                        raise TypeError  
                oldblock = block
                i2 += 1
                try:
                    line = self.parseline(line)
                except OverhangError :          #   overhang error 
                    n = name                    #   do not include the enzyme
                    if not bl[2]:
                        print 'Anyway, %s is not commercially available.\n' %n
                    else:
                        print 'Unfortunately, %s is commercially available.\n'%n

                    continue 
                #Hyphens can't be used as a Python name, nor as a
                #group name in a regular expression.
                name = name.replace("-","_")
                if name in enzymedict:
                    #
                    #   deal with TaqII and its two sites.
                    #
                    print '\nWARNING :',
                    print name, 'has two different sites.\n'
                    other = line[0].replace("-","_")
                    dna = DNA(line[1])
                    sense1 = regex(dna.tostring())
                    antisense1 = regex(Antiparallel(dna))
                    dna = DNA(enzymedict[other][0])
                    sense2 = regex(dna.tostring())
                    antisense2 = regex(Antiparallel(dna))
                    sense = '(?P<'+other+'>'+sense1+'|'+sense2+')'
                    antisense = '(?P<'+other+'_as>'+antisense1+'|'+antisense2 + ')'
                    reg = sense + '|' + antisense 
                    line[1] = line[1] + '|' + enzymedict[other][0]
                    line[-1] = reg
                #
                #   the data to produce the enzyme class are then stored in
                #   enzymedict.
                #
                enzymedict[name] = line[1:] #element zero was the name
        except IndexError:
            pass
        for i in supplier:
            #
            #   construction of the list of suppliers.
            #
            t = i.strip().split(' ', 1)
            suppliersdict[t[0]] = (t[1], [])
        return

    
