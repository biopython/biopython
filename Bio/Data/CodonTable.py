# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Codon tables based on those from the NCBI.

These tables are based on parsing the NCBI file:
ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt

Last updated for Version 3.9
"""

from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData

unambiguous_dna_by_name = {}
unambiguous_dna_by_id = {}
unambiguous_rna_by_name = {}
unambiguous_rna_by_id = {}
generic_by_name = {} # unambiguous DNA or RNA
generic_by_id = {} # unambiguous DNA or RNA

ambiguous_dna_by_name = {}
ambiguous_dna_by_id = {}
ambiguous_rna_by_name = {}
ambiguous_rna_by_id = {}
ambiguous_generic_by_name = {} # ambiguous DNA or RNA
ambiguous_generic_by_id = {} # ambiguous DNA or RNA 

# standard IUPAC unambiguous codons
standard_dna_table = None
standard_rna_table = None

# In the future, the back_table could return a statistically
# appropriate distribution of codons, so do not cache the results of
# back_table lookups!

class TranslationError(Exception):
    pass

class CodonTable(object):
    nucleotide_alphabet = Alphabet.generic_nucleotide
    protein_alphabet = Alphabet.generic_protein
    
    forward_table = {}    # only includes codons which actually code
    back_table = {}       # for back translations
    start_codons = []
    stop_codons = []
    # Not always called from derived classes!
    def __init__(self, nucleotide_alphabet = nucleotide_alphabet,
                 protein_alphabet = protein_alphabet,
                 forward_table = forward_table, back_table = back_table,
                 start_codons = start_codons, stop_codons = stop_codons):
        self.nucleotide_alphabet = nucleotide_alphabet
        self.protein_alphabet = protein_alphabet
        self.forward_table = forward_table
        self.back_table = back_table
        self.start_codons = start_codons
        self.stop_codons = stop_codons

    def __str__(self):
        """Returns a simple text representation of the codon table

        e.g.
        >>> import Bio.Data.CodonTable
        >>> print Bio.Data.CodonTable.standard_dna_table
        >>> print Bio.Data.CodonTable.generic_by_id[1]
        """

        if self.id:
            answer = "Table %i" % self.id
        else:
            answer = "Table ID unknown"
        if self.names:
            answer += " " + ", ".join(filter(None, self.names))

        #Use the main four letters (and the conventional ordering)
        #even for ambiguous tables
        letters = self.nucleotide_alphabet.letters
        if isinstance(self.nucleotide_alphabet, Alphabet.DNAAlphabet) \
        or (letters is not None and "T" in letters):
            letters = "TCAG"
        else:
            #Should be either RNA or generic nucleotides,
            #e.g. Bio.Data.CodonTable.generic_by_id[1]
            letters = "UCAG"

        #Build the table...
        answer=answer + "\n\n  |" + "|".join( \
            ["  %s      " % c2 for c2 in letters] \
            ) + "|"
        answer=answer + "\n--+" \
               + "+".join(["---------" for c2 in letters]) + "+--"
        for c1 in letters:
            for c3 in letters:
                line = c1 + " |"
                for c2 in letters:
                    codon = c1+c2+c3
                    line = line + " %s" % codon
                    if codon in self.stop_codons:
                        line = line + " Stop|"
                    else:
                        try:
                            amino = self.forward_table[codon]
                        except KeyError:
                            amino = "?"
                        except TranslationError:
                            amino = "?"
                        if codon in self.start_codons:
                            line = line + " %s(s)|" % amino
                        else:
                            line = line + " %s   |" % amino
                line = line + " " + c3
                answer = answer + "\n"+ line 
            answer=answer + "\n--+" \
                  + "+".join(["---------" for c2 in letters]) + "+--"
        return answer
            
def make_back_table(table, default_stop_codon):
    #  ONLY RETURNS A SINGLE CODON
    # Do the sort so changes in the hash implementation won't affect
    # the result when one amino acid is coded by more than one codon.
    back_table = {}
    for key in sorted(table):
        back_table[table[key]] = key
    back_table[None] = default_stop_codon
    return back_table


class NCBICodonTable(CodonTable):
    nucleotide_alphabet = Alphabet.generic_nucleotide
    protein_alphabet = IUPAC.protein
    
    def __init__(self, id, names, table, start_codons, stop_codons):
        self.id = id
        self.names = names
        self.forward_table = table
        self.back_table = make_back_table(table, stop_codons[0])
        self.start_codons = start_codons
        self.stop_codons = stop_codons


class NCBICodonTableDNA(NCBICodonTable):
    nucleotide_alphabet = IUPAC.unambiguous_dna

class NCBICodonTableRNA(NCBICodonTable):
    nucleotide_alphabet = IUPAC.unambiguous_rna


#########  Deal with ambiguous forward translations

class AmbiguousCodonTable(CodonTable):
    def __init__(self, codon_table,
                 ambiguous_nucleotide_alphabet,
                 ambiguous_nucleotide_values,
                 ambiguous_protein_alphabet,
                 ambiguous_protein_values):
        CodonTable.__init__(self,
                            ambiguous_nucleotide_alphabet,
                            ambiguous_protein_alphabet,
                            AmbiguousForwardTable(codon_table.forward_table,
                                                  ambiguous_nucleotide_values,
                                                  ambiguous_protein_values),
                            codon_table.back_table,

                            # These two are WRONG!  I need to get the
                            # list of ambiguous codons which code for
                            # the stop codons  XXX
                            list_ambiguous_codons(codon_table.start_codons, ambiguous_nucleotide_values),
                            list_ambiguous_codons(codon_table.stop_codons, ambiguous_nucleotide_values)
                            )
        self._codon_table = codon_table

    # Be sneaky and forward attribute lookups to the original table.
    # This lets us get the names, if the original table is an NCBI
    # table.
    def __getattr__(self, name):
        return getattr(self._codon_table, name)

def list_possible_proteins(codon, forward_table, ambiguous_nucleotide_values):
        c1, c2, c3 = codon
        x1 = ambiguous_nucleotide_values[c1]
        x2 = ambiguous_nucleotide_values[c2]
        x3 = ambiguous_nucleotide_values[c3]
        possible = {}
        stops = []
        for y1 in x1:
            for y2 in x2:
                for y3 in x3:
                    try:
                        possible[forward_table[y1+y2+y3]] = 1
                    except KeyError:
                        # If tripping over a stop codon
                        stops.append(y1+y2+y3)
        if stops:
            if possible:
                raise TranslationError("ambiguous codon '%s' codes " % codon \
                                       + "for both proteins and stop codons")
            # This is a true stop codon - tell the caller about it
            raise KeyError(codon)
        return possible.keys()

def list_ambiguous_codons(codons, ambiguous_nucleotide_values):
    """Extends a codon list to include all possible ambigous codons.

    e.g. ['TAG', 'TAA'] -> ['TAG', 'TAA', 'TAR']
         ['UAG', 'UGA'] -> ['UAG', 'UGA', 'URA']

    Note that ['TAG', 'TGA'] -> ['TAG', 'TGA'], this does not add 'TRR'.
    Thus only two more codons are added in the following:

    e.g. ['TGA', 'TAA', 'TAG'] -> ['TGA', 'TAA', 'TAG', 'TRA', 'TAR']

    Returns a new (longer) list of codon strings.
    """

    #Note ambiguous_nucleotide_values['R'] = 'AG' (etc)
    #This will generate things like 'TRR' from ['TAG', 'TGA'], which
    #we don't want to include:
    c1_list = sorted(letter for (letter, meanings) \
               in ambiguous_nucleotide_values.iteritems() \
               if set([codon[0] for codon in codons]).issuperset(set(meanings)))
    c2_list = sorted(letter for (letter, meanings) \
               in ambiguous_nucleotide_values.iteritems() \
               if set([codon[1] for codon in codons]).issuperset(set(meanings)))
    c3_list = sorted(letter for (letter, meanings) \
               in ambiguous_nucleotide_values.iteritems() \
               if set([codon[2] for codon in codons]).issuperset(set(meanings)))
    #candidates is a list (not a set) to preserve the iteration order
    candidates = []
    for c1 in c1_list:
        for c2 in c2_list:
            for c3 in c3_list:
                codon = c1+c2+c3
                if codon not in candidates and codon not in codons:
                    candidates.append(codon)
    answer = codons[:] #copy
    #print "Have %i new candidates" % len(candidates)
    for ambig_codon in candidates:
        wanted = True
        #e.g. 'TRR' -> 'TAA', 'TAG', 'TGA', 'TGG'
        for codon in [c1+c2+c3 \
                      for c1 in ambiguous_nucleotide_values[ambig_codon[0]] \
                      for c2 in ambiguous_nucleotide_values[ambig_codon[1]] \
                      for c3 in ambiguous_nucleotide_values[ambig_codon[2]]]:
            if codon not in codons:
                #This ambiguous codon can code for a non-stop, exclude it!
                wanted=False
                #print "Rejecting %s" % ambig_codon
                continue
        if wanted:
            answer.append(ambig_codon)
    return answer

assert list_ambiguous_codons(['TGA', 'TAA'],IUPACData.ambiguous_dna_values) == ['TGA', 'TAA', 'TRA']
assert list_ambiguous_codons(['TAG', 'TGA'],IUPACData.ambiguous_dna_values) == ['TAG', 'TGA']
assert list_ambiguous_codons(['TAG', 'TAA'],IUPACData.ambiguous_dna_values) == ['TAG', 'TAA', 'TAR']
assert list_ambiguous_codons(['UAG', 'UAA'],IUPACData.ambiguous_rna_values) == ['UAG', 'UAA', 'UAR']
assert list_ambiguous_codons(['TGA', 'TAA', 'TAG'],IUPACData.ambiguous_dna_values) == ['TGA', 'TAA', 'TAG', 'TAR', 'TRA']

# Forward translation is "onto", that is, any given codon always maps
# to the same protein, or it doesn't map at all.  Thus, I can build
# off of an existing table to produce the ambiguous mappings.
#
# This handles the general case.  Perhaps it's overkill?
#  >>> t = CodonTable.ambiguous_dna_by_id[1]
#  >>> t.forward_table["AAT"]
#  'N'
#  >>> t.forward_table["GAT"]
#  'D'
#  >>> t.forward_table["RAT"]
#  'B'
#  >>> t.forward_table["YTA"]
#  'L'

class AmbiguousForwardTable(object):
    def __init__(self, forward_table, ambiguous_nucleotide, ambiguous_protein):
        self.forward_table = forward_table

        self.ambiguous_nucleotide = ambiguous_nucleotide
        self.ambiguous_protein = ambiguous_protein

        inverted = {}
        for name, val in ambiguous_protein.iteritems():
            for c in val:
                x = inverted.get(c, {})
                x[name] = 1
                inverted[c] = x
        for name, val in inverted.iteritems():
            inverted[name] = val.keys()
        self._inverted = inverted
        
        self._cache = {}

    def get(self, codon, failobj = None):
        try:
            return self.__getitem__(codon)
        except KeyError:
            return failobj
        
    def __getitem__(self, codon):
        try:
            x = self._cache[codon]
        except KeyError:
            pass
        else:
            if x is TranslationError:
                raise TranslationError(codon)   # no unique translation
            if x is KeyError:
                raise KeyError(codon)  # it's a stop codon
            return x
        try:
            x = self.forward_table[codon]
            self._cache[codon] = x
            return x
        except KeyError:
            pass

        # XXX Need to make part of this into a method which returns
        # a list of all possible encodings for a codon!
        try:
            possible = list_possible_proteins(codon,
                                              self.forward_table,
                                              self.ambiguous_nucleotide)
        except KeyError:
            self._cache[codon] = KeyError
            raise KeyError(codon)  # stop codon
        except TranslationError:
            self._cache[codon] = TranslationError
            raise TranslationError(codon)  # does not code
        assert len(possible) > 0, "unambiguous codons must code"

        # Hah!  Only one possible protein, so use it
        if len(possible) == 1:
            self._cache[codon] = possible[0]
            return possible[0]

        # See if there's an ambiguous protein encoding for the multiples.
        # Find residues which exist in every coding set.
        ambiguous_possible = {}
        for amino in possible:
            for term in self._inverted[amino]:
                ambiguous_possible[term] = ambiguous_possible.get(term, 0) + 1

        n = len(possible)
        possible = []
        for amino, val in ambiguous_possible.iteritems():
            if val == n:
                possible.append(amino)

        # No amino acid encoding for the results
        if len(possible) == 0:
            self._cache[codon] = TranslationError
            raise TranslationError(codon)   # no valid translation

        # All of these are valid, so choose one
        # To be unique, sort by smallet ambiguity then alphabetically
        # Can get this if "X" encodes for everything.
        #def _sort(x, y, table = self.ambiguous_protein):
        #    a = cmp(len(table[x]), len(table[y]))
        #    if a == 0:
        #        return cmp(x, y)
        #    return a

        #Sort by key is 2.x and 3.x compatible
        possible.sort(key=lambda x:(len(self.ambiguous_protein[x]), x))
                          
        x = possible[0]
        self._cache[codon] = x
        return x


def register_ncbi_table(name, alt_name, id,
                        table, start_codons, stop_codons):
    """Turns codon table data into objects, and stores them in the dictionaries (PRIVATE)."""
    #In most cases names are divided by "; ", however there is also
    #'Bacterial and Plant Plastid' (which used to be just 'Bacterial')
    names = [x.strip() for x in name.replace(" and ","; ").split("; ")]
    
    dna = NCBICodonTableDNA(id, names + [alt_name], table, start_codons,
                            stop_codons)

    ambig_dna = AmbiguousCodonTable(dna,
                                    IUPAC.ambiguous_dna,
                                    IUPACData.ambiguous_dna_values,
                                    IUPAC.extended_protein,
                                    IUPACData.extended_protein_values)
    
    # replace all T's with U's for the RNA tables
    rna_table = {}
    generic_table = {}
    for codon, val in table.iteritems():
        generic_table[codon] = val
        codon = codon.replace("T", "U")
        generic_table[codon] = val
        rna_table[codon] = val
    rna_start_codons = []
    generic_start_codons = []
    for codon in start_codons:
        generic_start_codons.append(codon)
        codon = codon.replace("T", "U")
        generic_start_codons.append(codon)
        rna_start_codons.append(codon)
    rna_stop_codons = []
    generic_stop_codons = []
    for codon in stop_codons:
        generic_stop_codons.append(codon)
        codon = codon.replace("T", "U")
        generic_stop_codons.append(codon)
        rna_stop_codons.append(codon)
    
    generic = NCBICodonTable(id, names + [alt_name], generic_table,
                             generic_start_codons, generic_stop_codons)

    #The following isn't very elegant, but seems to work nicely.
    _merged_values = dict(IUPACData.ambiguous_rna_values.iteritems())
    _merged_values["T"] = "U"
    ambig_generic = AmbiguousCodonTable(generic,
                                        Alphabet.NucleotideAlphabet(),
                                        _merged_values,
                                        IUPAC.extended_protein,
                                        IUPACData.extended_protein_values)

    rna = NCBICodonTableRNA(id, names + [alt_name], rna_table,
                            rna_start_codons, rna_stop_codons)

    ambig_rna = AmbiguousCodonTable(rna,
                                    IUPAC.ambiguous_rna,
                                    IUPACData.ambiguous_rna_values,
                                    IUPAC.extended_protein,
                                    IUPACData.extended_protein_values)

    if id == 1:
        global standard_dna_table, standard_rna_table
        standard_dna_table = dna
        standard_rna_table = rna

    unambiguous_dna_by_id[id] = dna
    unambiguous_rna_by_id[id] = rna
    generic_by_id[id] = generic
    ambiguous_dna_by_id[id] = ambig_dna
    ambiguous_rna_by_id[id] = ambig_rna
    ambiguous_generic_by_id[id] = ambig_generic

    if alt_name is not None:
        names.append(alt_name)

    for name in names:
        unambiguous_dna_by_name[name] = dna
        unambiguous_rna_by_name[name] = rna
        generic_by_name[name] = generic
        ambiguous_dna_by_name[name] = ambig_dna
        ambiguous_rna_by_name[name] = ambig_rna
        ambiguous_generic_by_name[name] = ambig_generic


### These tables created from the data file
###  ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
### using the following:
##import re
##for line in open("gc.prt").readlines():
##    if line[:2] == " {":
##        names = []
##        id = None
##        aa = None
##        start = None
##        bases = []
##    elif line[:6] == "  name":
##        names.append(re.search('"([^"]*)"', line).group(1))
##    elif line[:8] == "    name":
##        names.append(re.search('"(.*)$', line).group(1))
##    elif line == ' Mitochondrial; Mycoplasma; Spiroplasma" ,\n':
##        names[-1] = names[-1] + " Mitochondrial; Mycoplasma; Spiroplasma"
##    elif line[:4] == "  id":
##        id = int(re.search('(\d+)', line).group(1))
##    elif line[:10] == "  ncbieaa ":
##        aa = line[12:12+64]
##    elif line[:10] == "  sncbieaa":
##        start = line[12:12+64]
##    elif line[:9] == "  -- Base":
##        bases.append(line[12:12+64])
##    elif line[:2] == " }":
##        assert names != [] and id is not None and aa is not None
##        assert start is not None and bases != []
##        if len(names) == 1:
##            names.append(None)
##        print "register_ncbi_table(name = %s," % repr(names[0])
##        print "                    alt_name = %s, id = %d," % \
##              (repr(names[1]), id)
##        print "                    table = {"
##        s = "    "
##        for i in range(64):
##            if aa[i] != "*":
##                t = " '%s%s%s': '%s'," % (bases[0][i], bases[1][i],
##                                          bases[2][i], aa[i])
##                if len(s) + len(t) > 75:
##                    print s
##                    s = "    " + t
##                else:
##                    s = s + t
##        print s, "},"

##        s = "                    stop_codons = ["
##        for i in range(64):
##            if aa[i] == "*":
##                t = " '%s%s%s'," % (bases[0][i], bases[1][i], bases[2][i])
##                if len(s) + len(t) > 75:
##                    print s
##                    s = "                                    " + t
##                else:
##                    s = s + t
##        print s, "],"

##        s = "                    start_codons = ["
##        for i in range(64):
##            if start[i] == "M":
##                t = " '%s%s%s'," % (bases[0][i], bases[1][i], bases[2][i])
##                if len(s) + len(t) > 75:
##                    print s
##                    s = "                                    " + t
##                else:
##                    s = s + t
##        print s, "]"
##        print "                    )"
##    elif line[:2] == "--" or line == "\n" or line == "}\n" or \
##         line == 'Genetic-code-table ::= {\n':
##        pass
##    else:
##        raise Exception("Unparsed: " + repr(line))

register_ncbi_table(name = 'Standard',
                    alt_name = 'SGC0', id = 1,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
     'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
     'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
     'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
     'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
     'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
     'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
     'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
     'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
     'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
     'GGG': 'G', },
                    stop_codons = [ 'TAA', 'TAG', 'TGA', ],
                    start_codons = [ 'TTG', 'CTG', 'ATG', ]
                    )
register_ncbi_table(name = 'Vertebrate Mitochondrial',
                    alt_name = 'SGC1', id = 2,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W', 'CTT': 'L',
     'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P',
     'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
     'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
     'ATT': 'I', 'ATC': 'I', 'ATA': 'M', 'ATG': 'M', 'ACT': 'T',
     'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
     'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'GTT': 'V',
     'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A',
     'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E',
     'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', },
                    stop_codons = [ 'TAA', 'TAG', 'AGA', 'AGG', ],
                    start_codons = [ 'ATT', 'ATC', 'ATA', 'ATG', 'GTG', ]
                    )
register_ncbi_table(name = 'Yeast Mitochondrial',
                    alt_name = 'SGC2', id = 3,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W', 'CTT': 'T',
     'CTC': 'T', 'CTA': 'T', 'CTG': 'T', 'CCT': 'P', 'CCC': 'P',
     'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
     'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
     'ATT': 'I', 'ATC': 'I', 'ATA': 'M', 'ATG': 'M', 'ACT': 'T',
     'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
     'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R',
     'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
     'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D',
     'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G',
     'GGA': 'G', 'GGG': 'G', },
                    stop_codons = [ 'TAA', 'TAG', ],
                    start_codons = [ 'ATA', 'ATG', ]
                    )
register_ncbi_table(name = 'Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma',
                    alt_name = 'SGC3', id = 4,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W', 'CTT': 'L',
     'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P',
     'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
     'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
     'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T',
     'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
     'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R',
     'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
     'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D',
     'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G',
     'GGA': 'G', 'GGG': 'G', },
                    stop_codons = [ 'TAA', 'TAG', ],
                    start_codons = [ 'TTA', 'TTG', 'CTG', 'ATT', 'ATC',
                                     'ATA', 'ATG', 'GTG', ]
                    )
register_ncbi_table(name = 'Invertebrate Mitochondrial',
                    alt_name = 'SGC4', id = 5,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W', 'CTT': 'L',
     'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P',
     'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
     'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
     'ATT': 'I', 'ATC': 'I', 'ATA': 'M', 'ATG': 'M', 'ACT': 'T',
     'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
     'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'S',
     'AGG': 'S', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
     'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D',
     'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G',
     'GGA': 'G', 'GGG': 'G', },
                    stop_codons = [ 'TAA', 'TAG', ],
                    start_codons = [ 'TTG', 'ATT', 'ATC', 'ATA', 'ATG',
                                     'GTG', ]
                    )
register_ncbi_table(name = 'Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear',
                    alt_name = 'SGC5', id = 6,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TAA': 'Q', 'TAG': 'Q', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W',
     'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P',
     'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H',
     'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
     'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
     'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N',
     'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S',
     'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
     'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
     'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G',
     'GGC': 'G', 'GGA': 'G', 'GGG': 'G', },
                    stop_codons = [ 'TGA', ],
                    start_codons = [ 'ATG', ]
                    )
register_ncbi_table(name = 'Echinoderm Mitochondrial; Flatworm Mitochondrial',
                    alt_name = 'SGC8', id = 9,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W', 'CTT': 'L',
     'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P',
     'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
     'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
     'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T',
     'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
     'AAA': 'N', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'S',
     'AGG': 'S', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
     'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D',
     'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G',
     'GGA': 'G', 'GGG': 'G', },
                    stop_codons = [ 'TAA', 'TAG', ],
                    start_codons = [ 'ATG', 'GTG', ]
                    )
register_ncbi_table(name = 'Euplotid Nuclear',
                    alt_name = 'SGC9', id = 10,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TGT': 'C', 'TGC': 'C', 'TGA': 'C', 'TGG': 'W', 'CTT': 'L',
     'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P',
     'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
     'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
     'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T',
     'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
     'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R',
     'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
     'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D',
     'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G',
     'GGA': 'G', 'GGG': 'G', },
                    stop_codons = [ 'TAA', 'TAG', ],
                    start_codons = [ 'ATG', ]
                    )
register_ncbi_table(name = 'Bacterial and Plant Plastid',
                    alt_name = None, id = 11,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
     'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
     'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
     'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
     'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
     'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
     'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
     'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
     'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
     'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
     'GGG': 'G', },
                    stop_codons = [ 'TAA', 'TAG', 'TGA', ],
                    start_codons = [ 'TTG', 'CTG', 'ATT', 'ATC', 'ATA',
                                     'ATG', 'GTG', ]
                    )
register_ncbi_table(name = 'Alternative Yeast Nuclear',
                    alt_name = None, id = 12,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
     'CTA': 'L', 'CTG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
     'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
     'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
     'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
     'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
     'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
     'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
     'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
     'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
     'GGG': 'G', },
                    stop_codons = [ 'TAA', 'TAG', 'TGA', ],
                    start_codons = [ 'CTG', 'ATG', ]
                    )
register_ncbi_table(name = 'Ascidian Mitochondrial',
                    alt_name = None, id = 13,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W', 'CTT': 'L',
     'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P',
     'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
     'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
     'ATT': 'I', 'ATC': 'I', 'ATA': 'M', 'ATG': 'M', 'ACT': 'T',
     'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
     'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'G',
     'AGG': 'G', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
     'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D',
     'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G',
     'GGA': 'G', 'GGG': 'G', },
                    stop_codons = [ 'TAA', 'TAG', ],
                    start_codons = [ 'TTG', 'ATA', 'ATG', 'GTG', ]
                    )
register_ncbi_table(name = 'Alternative Flatworm Mitochondrial',
                    alt_name = None, id = 14,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TAA': 'Y', 'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W',
     'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P',
     'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H',
     'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
     'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
     'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N',
     'AAC': 'N', 'AAA': 'N', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S',
     'AGA': 'S', 'AGG': 'S', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
     'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
     'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G',
     'GGC': 'G', 'GGA': 'G', 'GGG': 'G', },
                    stop_codons = [ 'TAG', ],
                    start_codons = [ 'ATG', ]
                    )
register_ncbi_table(name = 'Blepharisma Macronuclear',
                    alt_name = None, id = 15,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TAG': 'Q', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L',
     'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P',
     'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
     'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
     'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T',
     'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
     'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R',
     'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
     'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D',
     'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G',
     'GGA': 'G', 'GGG': 'G', },
                    stop_codons = [ 'TAA', 'TGA', ],
                    start_codons = [ 'ATG', ]
                    )
register_ncbi_table(name = 'Chlorophycean Mitochondrial',
                    alt_name = None, id = 16,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TAG': 'L', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L',
     'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P',
     'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
     'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
     'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T',
     'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
     'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R',
     'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
     'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D',
     'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G',
     'GGA': 'G', 'GGG': 'G', },
                    stop_codons = [ 'TAA', 'TGA', ],
                    start_codons = [ 'ATG', ]
                    )
register_ncbi_table(name = 'Trematode Mitochondrial',
                    alt_name = None, id = 21,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W', 'CTT': 'L',
     'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P',
     'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
     'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
     'ATT': 'I', 'ATC': 'I', 'ATA': 'M', 'ATG': 'M', 'ACT': 'T',
     'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
     'AAA': 'N', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'S',
     'AGG': 'S', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
     'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D',
     'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G',
     'GGA': 'G', 'GGG': 'G', },
                    stop_codons = [ 'TAA', 'TAG', ],
                    start_codons = [ 'ATG', 'GTG', ]
                    )
register_ncbi_table(name = 'Scenedesmus obliquus Mitochondrial',
                    alt_name = None, id = 22,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAG': 'L',
     'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
     'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
     'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
     'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
     'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
     'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
     'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
     'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
     'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
     'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
     'GGG': 'G', },
                    stop_codons = [ 'TCA', 'TAA', 'TGA', ],
                    start_codons = [ 'ATG', ]
                    )
register_ncbi_table(name = 'Thraustochytrium Mitochondrial',
                    alt_name = None, id = 23,
                    table = {
     'TTT': 'F', 'TTC': 'F', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S',
     'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C',
     'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L',
     'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
     'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R',
     'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I',
     'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T',
     'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
     'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V',
     'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A',
     'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E',
     'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', },
                    stop_codons = [ 'TTA', 'TAA', 'TAG', 'TGA', ],
                    start_codons = [ 'ATT', 'ATG', 'GTG', ]
                    )



#Basic sanity test,
for key, val in generic_by_name.iteritems():
    assert key in ambiguous_generic_by_name[key].names
for key, val in generic_by_id.iteritems():
    assert ambiguous_generic_by_id[key].id == key
del key, val

for n in ambiguous_generic_by_id:
    assert ambiguous_rna_by_id[n].forward_table["GUU"] == "V"
    assert ambiguous_rna_by_id[n].forward_table["GUN"] == "V"
    if n != 23 :
        #For table 23, UUN = F, L or stop.
        assert ambiguous_rna_by_id[n].forward_table["UUN"] == "X" #F or L
    #R = A or G, so URR = UAA or UGA / TRA = TAA or TGA = stop codons
    if "UAA" in unambiguous_rna_by_id[n].stop_codons \
    and "UGA" in unambiguous_rna_by_id[n].stop_codons:
        try:
            print ambiguous_dna_by_id[n].forward_table["TRA"]
            assert False, "Should be a stop only"
        except KeyError:
            pass
        assert "URA" in ambiguous_generic_by_id[n].stop_codons
        assert "URA" in ambiguous_rna_by_id[n].stop_codons
        assert "TRA" in ambiguous_generic_by_id[n].stop_codons
        assert "TRA" in ambiguous_dna_by_id[n].stop_codons
del n
assert ambiguous_generic_by_id[1] == ambiguous_generic_by_name["Standard"]
assert ambiguous_generic_by_id[4] == ambiguous_generic_by_name["SGC3"]
assert ambiguous_generic_by_id[11] == ambiguous_generic_by_name["Bacterial"]
assert ambiguous_generic_by_id[11] == ambiguous_generic_by_name["Plant Plastid"]
assert ambiguous_generic_by_id[15] == ambiguous_generic_by_name['Blepharisma Macronuclear']
assert generic_by_id[1] == generic_by_name["Standard"]
assert generic_by_id[4] == generic_by_name["SGC3"]
assert generic_by_id[11] == generic_by_name["Bacterial"]
assert generic_by_id[11] == generic_by_name["Plant Plastid"]
assert generic_by_id[15] == generic_by_name['Blepharisma Macronuclear']
