# Copyright 2007-2009 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package..

"""Bio.SeqIO support for the "genbank" and "embl" file formats.

You are expected to use this module via the Bio.SeqIO functions.
Note that internally this module calls Bio.GenBank to do the actual
parsing of both GenBank and EMBL files.

See also:

International Nucleotide Sequence Database Collaboration
http://www.insdc.org/
 
GenBank
http://www.ncbi.nlm.nih.gov/Genbank/

EMBL Nucleotide Sequence Database
http://www.ebi.ac.uk/embl/

DDBJ (DNA Data Bank of Japan)
http://www.ddbj.nig.ac.jp/
"""

from Bio.Seq import UnknownSeq
from Bio.GenBank.Scanner import GenBankScanner, EmblScanner
from Bio import Alphabet
from Interfaces import SequentialSequenceWriter

# NOTE
# ====
# The "brains" for parsing GenBank and EMBL files (and any
# other flat file variants from the INSDC in future) is in
# Bio.GenBank.Scanner (plus the _FeatureConsumer in Bio.GenBank)

def GenBankIterator(handle) :
    """Breaks up a Genbank file into SeqRecord objects.

    Every section from the LOCUS line to the terminating // becomes
    a single SeqRecord with associated annotation and features.
    
    Note that for genomes or chromosomes, there is typically only
    one record."""
    #This calls a generator function:
    return GenBankScanner(debug=0).parse_records(handle)

def EmblIterator(handle) :
    """Breaks up an EMBL file into SeqRecord objects.

    Every section from the LOCUS line to the terminating // becomes
    a single SeqRecord with associated annotation and features.
    
    Note that for genomes or chromosomes, there is typically only
    one record."""
    #This calls a generator function:
    return EmblScanner(debug=0).parse_records(handle)

def GenBankCdsFeatureIterator(handle, alphabet=Alphabet.generic_protein) :
    """Breaks up a Genbank file into SeqRecord objects for each CDS feature.

    Every section from the LOCUS line to the terminating // can contain
    many CDS features.  These are returned as with the stated amino acid
    translation sequence (if given).
    """
    #This calls a generator function:
    return GenBankScanner(debug=0).parse_cds_features(handle, alphabet)
    
def EmblCdsFeatureIterator(handle, alphabet=Alphabet.generic_protein) :
    """Breaks up a EMBL file into SeqRecord objects for each CDS feature.

    Every section from the LOCUS line to the terminating // can contain
    many CDS features.  These are returned as with the stated amino acid
    translation sequence (if given).
    """
    #This calls a generator function:
    return EmblScanner(debug=0).parse_cds_features(handle, alphabet)

class GenBankWriter(SequentialSequenceWriter) :
    HEADER_WIDTH = 12
    MAX_WIDTH = 80
    
    def _write_single_line(self, tag, text) :
        "Used in the the 'header' of each GenBank record."""
        assert len(tag) < self.HEADER_WIDTH
        assert len(text) < self.MAX_WIDTH - self.HEADER_WIDTH, \
               "Annotation %s too long for %s line" % (repr(text), tag)
        self.handle.write("%s%s\n" % (tag.ljust(self.HEADER_WIDTH),
                                      text.replace("\n"," ")))

    def _write_multi_line(self, tag, text) :
        "Used in the the 'header' of each GenBank record."""
        #TODO - Do the line spliting while preserving white space?
        max_len = self.MAX_WIDTH - self.HEADER_WIDTH
        assert len(tag) < self.HEADER_WIDTH
        text = text.strip()
        if len(text) < max_len :
            self._write_single_line(tag, text)
            return

        words = text.split()
        assert max([len(w) for w in words]) < max_len, \
               "Your description cannot be broken into nice lines!"
        text = ""
        while words and len(text) + 1 + len(words[0]) < max_len :
            text += " " + words.pop(0)
            text = text.strip()
        assert len(text) < max_len
        self._write_single_line(tag, text)
        while words :
            text = ""
            while words and len(text) + 1 + len(words[0]) < max_len :
                text += " " + words.pop(0)
                text = text.strip()
            assert len(text) < max_len
            self._write_single_line("", text)
        assert not words

    def _write_the_first_line(self, record) :
        """Write the LOCUS line."""
        
        locus = record.name
        if not locus or locus == "<unknown name>" :
            locus = record.id
        if not locus or locus == "<unknown id>" :
            locus = self._get_annotation_str(record, "accession", just_first=True)
        if len(locus) > 16 :
            raise ValueError("Locus identifier %s is too long" % repr(locus))

        if len(record) > 99999999999 :
            #Currently GenBank only officially support up to 350000, but
            #the length field can take eleven digits
            raise ValueError("Sequence too long!")

        #Get the base alphabet (underneath any Gapped or StopCodon encoding)
        a = Alphabet._get_base_alphabet(record.seq.alphabet)
        if not isinstance(a, Alphabet.Alphabet) :
            raise TypeError("Invalid alphabet")
        elif isinstance(a, Alphabet.ProteinAlphabet) :
            units = "bp"
        elif isinstance(a, Alphabet.NucleotideAlphabet) :
            units = "aa"
        else :
            #Must be something like NucleotideAlphabet or
            #just the generic Alphabet (default for fasta files)
            raise ValueError("Need a Nucleotide or Protein alphabet")

        #Get the molecule type
        #TODO - record this explicitly in the parser?
        if isinstance(a, Alphabet.ProteinAlphabet) :
            mol_type = ""
        elif isinstance(a, Alphabet.DNAAlphabet) :
            mol_type = "DNA"
        elif isinstance(a, Alphabet.RNAAlphabet) :
            mol_type = "RNA"
        else :
            #Must be something like NucleotideAlphabet or
            #just the generic Alphabet (default for fasta files)
            raise ValueError("Need a DNA, RNA or Protein alphabet")
        
        try :
            division = record.annotations["data_file_division"]
        except KeyError :
            division = "UNK"
        if division not in ["PRI","ROD","MAM","VRT","INV","PLN","BCT",
                            "VRL","PHG","SYN","UNA","EST","PAT","STS",
                            "GSS","HTG","HTC","ENV"] :
            division = "UNK"
        
        assert len(units) == 2
        assert len(division) == 3
        #TODO - date
        #TODO - mol_type
        line = "LOCUS       %s %s %s    %s           %s 01-JAN-1980\n" \
                     % (locus.ljust(16),
                        str(len(record)).rjust(11),
                        units,
                        mol_type.ljust(6),
                        division)
        assert len(line) == 79+1, repr(line) #plus one for new line

        assert line[12:28].rstrip() == locus, \
               'LOCUS line does not contain the locus at the expected position:\n' + line
        assert line[28:29] == " "
        assert line[29:40].lstrip() == str(len(record)), \
               'LOCUS line does not contain the length at the expected position:\n' + line

        #Tests copied from Bio.GenBank.Scanner
        assert line[40:44] in [' bp ', ' aa '] , \
               'LOCUS line does not contain size units at expected position:\n' + line
        assert line[44:47] in ['   ', 'ss-', 'ds-', 'ms-'], \
               'LOCUS line does not have valid strand type (Single stranded, ...):\n' + line
        assert line[47:54].strip() == "" \
        or line[47:54].strip().find('DNA') != -1 \
        or line[47:54].strip().find('RNA') != -1, \
               'LOCUS line does not contain valid sequence type (DNA, RNA, ...):\n' + line
        assert line[54:55] == ' ', \
               'LOCUS line does not contain space at position 55:\n' + line
        assert line[55:63].strip() in ['','linear','circular'], \
               'LOCUS line does not contain valid entry (linear, circular, ...):\n' + line
        assert line[63:64] == ' ', \
               'LOCUS line does not contain space at position 64:\n' + line
        assert line[67:68] == ' ', \
               'LOCUS line does not contain space at position 68:\n' + line
        assert line[70:71] == '-', \
               'LOCUS line does not contain - at position 71 in date:\n' + line
        assert line[74:75] == '-', \
               'LOCUS line does not contain - at position 75 in date:\n' + line

        self.handle.write(line)

    def _get_annotation_str(self, record, key, default=".", just_first=False) :
        """Get an annotation dictionary entry (as a string).

        Some entries are lists, in which case if just_first=True the first entry
        is returned.  If just_first=False (default) this verifies there is only
        one entry before returning it."""
        try :
            answer = record.annotations[key]
        except KeyError :
            return default
        if isinstance(answer, list) :
            if not just_first : assert len(answer) == 1
            return str(answer[0])
        else :
            return str(answer)

    def _write_sequence(self, record):
        #Loosely based on code from Howard Salis
        #TODO - Force lower case?
        LETTERS_PER_LINE = 60
        SEQUENCE_INDENT = 9

        if isinstance(record.seq, UnknownSeq) :
            #We have already recorded the length, and there is no need
            #to record a long sequence of NNNNNNN...NNN or whatever.
            return

        data = self._get_seq_string(record) #Catches sequence being None
        seq_len = len(data)
        for line_number in range(0,seq_len,LETTERS_PER_LINE):
            self.handle.write(str(line_number+1).rjust(SEQUENCE_INDENT))
            for words in range(line_number,min(line_number+LETTERS_PER_LINE,seq_len),10):
                self.handle.write(" %s" % data[words:words+10])
            self.handle.write("\n")
        
    def write_record(self, record):
        """Write a single record to the output file."""
        handle = self.handle
        self._write_the_first_line(record)

        accession = self._get_annotation_str(record, "accession",
                                             record.id.split(".",1)[0],
                                             just_first=True)
        acc_with_version = accession
        if record.id.startswith(accession+".") :
            try :
                acc_with_version = "%s.%i" \
                                   % (accession, int(record.id.split(".",1)[1]))
            except ValueError :
                pass
        gi = self._get_annotation_str(record, "gi", just_first=True)

        descr = record.description
        if descr == "<unknown description>" : descr = "."
        self._write_multi_line("DEFINITION", descr)
        
        self._write_single_line("ACCESSION", accession)
        if gi != "." :
            self._write_single_line("VERSION", "%s  GI:%s" % (acc_with_version,gi))
        else :
            self._write_single_line("VERSION", "%s" % (acc_with_version))

        try :
            #List of strings
            keywords = "; ".join(record.annotations["keywords"])
        except KeyError :
            keywords = "."
        self._write_multi_line("KEYWORDS", keywords)

        self._write_multi_line("SOURCE", \
                                self._get_annotation_str(record, "source"))
        #The ORGANISM line MUST be a single line, as any continuation is the taxonomy
        org = self._get_annotation_str(record, "organism")
        if len(org) > self.MAX_WIDTH - self.HEADER_WIDTH :
            org = org[:self.MAX_WIDTH - self.HEADER_WIDTH-4]+"..."
        self._write_single_line("  ORGANISM", org)
        try :
            #List of strings
            taxonomy = "; ".join(record.annotations["taxonomy"])
        except KeyError :
            taxonomy = "."
        self._write_multi_line("", taxonomy)

        #TODO - References...
        handle.write("FEATURES             Location/Qualifiers\n")
        for feature in record.features :
            self._write_feature(feature) 
        handle.write("ORIGIN\n")
        self._write_sequence(record)
        handle.write("//\n")

    def _write_feature(self, feature):
        """Write a single SeqFeature object to features table.

        Not implemented yet, but this stub exists in the short term to
        facilitate working on writing GenBank files with a sub-class."""
        #TODO - Features...
        pass

if __name__ == "__main__" :
    print "Quick self test"
    import os
    from StringIO import StringIO

    def check_genbank_writer(records) :
        handle = StringIO()
        GenBankWriter(handle).write_file(records)
        handle.seek(0)

        records2 = list(GenBankIterator(handle))

        assert len(records) == len(records2)
        for r1, r2 in zip(records, records2) :
            #The SwissProt parser may leave \n in the description...
            assert r1.description.replace("\n", " ") == r2.description
            assert r1.id == r2.id
            assert r1.name == r2.name
            assert str(r1.seq) == str(r2.seq)
            for key in ["gi", "keywords", "source", "taxonomy"] :
                if key in r1.annotations :
                    assert r1.annotations[key] == r2.annotations[key], key
            for key in ["organism"] :
                if key in r1.annotations :
                    v1 = r1.annotations[key]
                    v2 = r2.annotations[key]
                    assert isinstance(v1, str) and isinstance(v2, str)
                    #SwissProt organism can be too long to record in GenBank format
                    assert v1 == v2 or \
                           (v2.endswith("...") and v1.startswith(v2[:-3])), key

    for filename in os.listdir("../../Tests/GenBank") :
        if not filename.endswith(".gbk") and not filename.endswith(".gb") :
            continue
        print filename
        
        handle = open("../../Tests/GenBank/%s" % filename)
        records = list(GenBankIterator(handle))
        handle.close()

        check_genbank_writer(records)

    for filename in os.listdir("../../Tests/EMBL") :
        if not filename.endswith(".embl") :
            continue
        print filename
        
        handle = open("../../Tests/EMBL/%s" % filename)
        records = list(EmblIterator(handle))
        handle.close()

        check_genbank_writer(records)

    from Bio import SeqIO
    for filename in os.listdir("../../Tests/SwissProt") :
        if not filename.startswith("sp") :
            continue
        print filename
        
        handle = open("../../Tests/SwissProt/%s" % filename)
        records = list(SeqIO.parse(handle,"swiss"))
        handle.close()

        check_genbank_writer(records)

