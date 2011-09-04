# Copyright 2007-2011 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package..

"""Bio.SeqIO support for the "genbank" and "embl" file formats.

You are expected to use this module via the Bio.SeqIO functions.
Note that internally this module calls Bio.GenBank to do the actual
parsing of GenBank, EMBL and IMGT files.

See also:

International Nucleotide Sequence Database Collaboration
http://www.insdc.org/
 
GenBank
http://www.ncbi.nlm.nih.gov/Genbank/

EMBL Nucleotide Sequence Database
http://www.ebi.ac.uk/embl/

DDBJ (DNA Data Bank of Japan)
http://www.ddbj.nig.ac.jp/

IMGT (use a variant of EMBL format with longer feature indents)
http://imgt.cines.fr/download/LIGM-DB/userman_doc.html
http://imgt.cines.fr/download/LIGM-DB/ftable_doc.html
http://www.ebi.ac.uk/imgt/hla/docs/manual.html

"""

from Bio.Seq import UnknownSeq
from Bio.GenBank.Scanner import GenBankScanner, EmblScanner, _ImgtScanner
from Bio import Alphabet
from Interfaces import SequentialSequenceWriter
from Bio import SeqFeature

from Bio._py3k import _is_int_or_long

# NOTE
# ====
# The "brains" for parsing GenBank, EMBL and IMGT files (and any
# other flat file variants from the INSDC in future) is in
# Bio.GenBank.Scanner (plus the _FeatureConsumer in Bio.GenBank)
# However, all the writing code is in this file.


def GenBankIterator(handle):
    """Breaks up a Genbank file into SeqRecord objects.

    Every section from the LOCUS line to the terminating // becomes
    a single SeqRecord with associated annotation and features.
    
    Note that for genomes or chromosomes, there is typically only
    one record."""
    #This calls a generator function:
    return GenBankScanner(debug=0).parse_records(handle)

def EmblIterator(handle):
    """Breaks up an EMBL file into SeqRecord objects.

    Every section from the LOCUS line to the terminating // becomes
    a single SeqRecord with associated annotation and features.
    
    Note that for genomes or chromosomes, there is typically only
    one record."""
    #This calls a generator function:
    return EmblScanner(debug=0).parse_records(handle)

def ImgtIterator(handle):
    """Breaks up an IMGT file into SeqRecord objects.

    Every section from the LOCUS line to the terminating // becomes
    a single SeqRecord with associated annotation and features.
    
    Note that for genomes or chromosomes, there is typically only
    one record."""
    #This calls a generator function:
    return _ImgtScanner(debug=0).parse_records(handle)

def GenBankCdsFeatureIterator(handle, alphabet=Alphabet.generic_protein):
    """Breaks up a Genbank file into SeqRecord objects for each CDS feature.

    Every section from the LOCUS line to the terminating // can contain
    many CDS features.  These are returned as with the stated amino acid
    translation sequence (if given).
    """
    #This calls a generator function:
    return GenBankScanner(debug=0).parse_cds_features(handle, alphabet)
    
def EmblCdsFeatureIterator(handle, alphabet=Alphabet.generic_protein):
    """Breaks up a EMBL file into SeqRecord objects for each CDS feature.

    Every section from the LOCUS line to the terminating // can contain
    many CDS features.  These are returned as with the stated amino acid
    translation sequence (if given).
    """
    #This calls a generator function:
    return EmblScanner(debug=0).parse_cds_features(handle, alphabet)

def _insdc_feature_position_string(pos, offset=0):
    """Build a GenBank/EMBL position string (PRIVATE).

    Use offset=1 to add one to convert a start position from python counting.
    """
    if isinstance(pos, SeqFeature.ExactPosition):
        return "%i" % (pos.position+offset)
    elif isinstance(pos, SeqFeature.WithinPosition):
        return "(%i.%i)" % (pos.position + offset,
                            pos.position + pos.extension + offset)
    elif isinstance(pos, SeqFeature.BetweenPosition):
        return "(%i^%i)" % (pos.position + offset,
                            pos.position + pos.extension + offset)
    elif isinstance(pos, SeqFeature.BeforePosition):
        return "<%i" % (pos.position + offset)
    elif isinstance(pos, SeqFeature.AfterPosition):
        return ">%i" % (pos.position + offset)
    elif isinstance(pos, SeqFeature.OneOfPosition):
        return "one-of(%s)" \
               % ",".join([_insdc_feature_position_string(p,offset) \
                           for p in pos.position_choices])
    elif isinstance(pos, SeqFeature.AbstractPosition):
        raise NotImplementedError("Please report this as a bug in Biopython.")
    else:
        raise ValueError("Expected a SeqFeature position object.")


def _insdc_location_string_ignoring_strand_and_subfeatures(location, rec_length):
    if location.ref:
        ref = "%s:" % location.ref
    else:
        ref = ""
    assert not location.ref_db
    if isinstance(location.start, SeqFeature.ExactPosition) \
    and isinstance(location.end, SeqFeature.ExactPosition) \
    and location.start.position == location.end.position:
        #Special case, for 12:12 return 12^13
        #(a zero length slice, meaning the point between two letters)
        if location.end.position == rec_length:
            #Very special case, for a between position at the end of a
            #sequence (used on some circular genomes, Bug 3098) we have
            #N:N so return N^1
            return "%s%i^1" % (ref, rec_length)
        else:
            return "%s%i^%i" % (ref, location.end.position,
                                location.end.position+1)
    if isinstance(location.start, SeqFeature.ExactPosition) \
    and isinstance(location.end, SeqFeature.ExactPosition) \
    and location.start.position+1 == location.end.position:
        #Special case, for 11:12 return 12 rather than 12..12
        #(a length one slice, meaning a single letter)
        return "%s%i" % (ref, location.end.position)
    elif isinstance(location.start, SeqFeature.UnknownPosition) \
    or isinstance(location.end, SeqFeature.UnknownPosition):
        #Special case for features from SwissProt/UniProt files
        if isinstance(location.start, SeqFeature.UnknownPosition) \
        and isinstance(location.end, SeqFeature.UnknownPosition):
            #import warnings
            #warnings.warn("Feature with unknown location")
            #return "?"
            raise ValueError("Feature with unknown location")
        elif isinstance(location.start, SeqFeature.UnknownPosition):
            #Treat the unknown start position as a BeforePosition
            return "%s<%i..%s" \
                % (ref,
                   location.nofuzzy_end,
                   _insdc_feature_position_string(location.end))
        else:
            #Treat the unknown end position as an AfterPosition
            return "%s%s..>%i" \
                % (ref,
                   _insdc_feature_position_string(location.start),
                   location.nofuzzy_start)
    else:
        #Typical case, e.g. 12..15 gets mapped to 11:15
        return ref \
               + _insdc_feature_position_string(location.start, +1) \
               + ".." + \
               _insdc_feature_position_string(location.end)

def _insdc_feature_location_string(feature, rec_length):
    """Build a GenBank/EMBL location string from a SeqFeature (PRIVATE).

    There is a choice of how to show joins on the reverse complement strand,
    GenBank used "complement(join(1,10),(20,100))" while EMBL used to use
    "join(complement(20,100),complement(1,10))" instead (but appears to have
    now adopted the GenBank convention). Notice that the order of the entries
    is reversed! This function therefore uses the first form. In this situation
    we expect the parent feature and the two children to all be marked as
    strand == -1, and in the order 0:10 then 19:100.

    Also need to consider dual-strand examples like these from the Arabidopsis
    thaliana chloroplast NC_000932: join(complement(69611..69724),139856..140650)
    gene ArthCp047, GeneID:844801 or its CDS (protein NP_051038.1 GI:7525057)
    which is further complicated by a splice:
    join(complement(69611..69724),139856..140087,140625..140650)

    For mixed this mixed strand feature, the parent SeqFeature should have
    no strand (either 0 or None) while the child features should have either
    strand +1 or -1 as appropriate, and be listed in the order given here.
    """

    if not feature.sub_features:
        #Non-recursive.
        #assert feature.location_operator == "", \
        #       "%s has no subfeatures but location_operator %s" \
        #       % (repr(feature), feature.location_operator)
        location = _insdc_location_string_ignoring_strand_and_subfeatures(feature.location, rec_length)
        if feature.strand == -1:
            location = "complement(%s)" % location
        return location
    # As noted above, treat reverse complement strand features carefully:
    if feature.strand == -1:
        for f in feature.sub_features:
            if f.strand != -1:
                raise ValueError("Inconsistent strands: %r for parent, %r for child" \
                                 % (feature.strand, f.strand))
        return "complement(%s(%s))" \
               % (feature.location_operator,
                  ",".join(_insdc_location_string_ignoring_strand_and_subfeatures(f.location, rec_length) \
                           for f in feature.sub_features))
    #if feature.strand == +1:
    #    for f in feature.sub_features:
    #        assert f.strand == +1
    #This covers typical forward strand features, and also an evil mixed strand:
    assert feature.location_operator != ""
    return  "%s(%s)" % (feature.location_operator,
                        ",".join([_insdc_feature_location_string(f, rec_length) \
                                  for f in feature.sub_features]))


class _InsdcWriter(SequentialSequenceWriter):
    """Base class for GenBank and EMBL writers (PRIVATE)."""
    MAX_WIDTH = 80
    QUALIFIER_INDENT = 21
    QUALIFIER_INDENT_STR = " "*QUALIFIER_INDENT
    QUALIFIER_INDENT_TMP = "     %s                " # 21 if %s is empty

    def _write_feature_qualifier(self, key, value=None, quote=None):
        if not value:
            self.handle.write("%s/%s\n" % (self.QUALIFIER_INDENT_STR, key))
            return
        #Quick hack with no line wrapping, may be useful for testing:
        #self.handle.write('%s/%s="%s"\n' % (self.QUALIFIER_INDENT_STR, key, value))
        if quote is None:
            #Try to mimic unwritten rules about when quotes can be left out:
            if _is_int_or_long(value):
                quote = False
            else:
                quote = True
        if quote:
            line = '%s/%s="%s"' % (self.QUALIFIER_INDENT_STR, key, value)
        else:
            line = '%s/%s=%s' % (self.QUALIFIER_INDENT_STR, key, value)
        if len(line) <= self.MAX_WIDTH:
            self.handle.write(line+"\n")
            return
        while line.lstrip():
            if len(line) <= self.MAX_WIDTH:
                self.handle.write(line+"\n")
                return
            #Insert line break...
            for index in range(min(len(line)-1, self.MAX_WIDTH),
                               self.QUALIFIER_INDENT+1,-1):
                if line[index] == " " : break
            if line[index] != " ":
                #No nice place to break...
                index = self.MAX_WIDTH
            assert index <= self.MAX_WIDTH
            self.handle.write(line[:index] + "\n")
            line = self.QUALIFIER_INDENT_STR + line[index:].lstrip()

    def _wrap_location(self, location):
        """Split a feature location into lines (break at commas)."""
        #TODO - Rewrite this not to recurse!
        length = self.MAX_WIDTH - self.QUALIFIER_INDENT
        if len(location) <= length:
            return location
        index = location[:length].rfind(",")
        if index == -1:
            #No good place to split (!)
            import warnings
            warnings.warn("Couldn't split location:\n%s" % location)
            return location
        return location[:index+1] + "\n" + \
               self.QUALIFIER_INDENT_STR + self._wrap_location(location[index+1:])

    def _write_feature(self, feature, record_length):
        """Write a single SeqFeature object to features table."""
        assert feature.type, feature
        location = _insdc_feature_location_string(feature, record_length)
        f_type = feature.type.replace(" ","_")
        line = (self.QUALIFIER_INDENT_TMP  % f_type)[:self.QUALIFIER_INDENT] \
               + self._wrap_location(location) + "\n"
        self.handle.write(line)
        #Now the qualifiers...
        for key, values in feature.qualifiers.iteritems():
            if isinstance(values, list) or isinstance(values, tuple):
                for value in values:
                    self._write_feature_qualifier(key, value)
            elif values:
                #String, int, etc
                self._write_feature_qualifier(key, values)
            else:
                #e.g. a /psuedo entry
                self._write_feature_qualifier(key)

    def _get_annotation_str(self, record, key, default=".", just_first=False):
        """Get an annotation dictionary entry (as a string).

        Some entries are lists, in which case if just_first=True the first entry
        is returned.  If just_first=False (default) this verifies there is only
        one entry before returning it."""
        try:
            answer = record.annotations[key]
        except KeyError:
            return default
        if isinstance(answer, list):
            if not just_first : assert len(answer) == 1
            return str(answer[0])
        else:
            return str(answer)

    def _split_multi_line(self, text, max_len):
        """Returns a list of strings.
        
        Any single words which are too long get returned as a whole line
        (e.g. URLs) without an exception or warning.
        """
        #TODO - Do the line spliting while preserving white space?
        text = text.strip()
        if len(text) <= max_len:
            return [text]

        words = text.split()
        text = ""
        while words and len(text) + 1 + len(words[0]) <= max_len:
            text += " " + words.pop(0)
            text = text.strip()
        #assert len(text) <= max_len
        answer = [text]
        while words:
            text = words.pop(0)
            while words and len(text) + 1 + len(words[0]) <= max_len:
                text += " " + words.pop(0)
                text = text.strip()
            #assert len(text) <= max_len
            answer.append(text)
        assert not words
        return answer

    def _split_contig(self, record, max_len):
        "Returns a list of strings, splits on commas."""
        #TODO - Merge this with _write_multi_line method?
        #It would need the addition of the comma splitting logic...
        #are there any other cases where that would be sensible?
        contig = record.annotations.get("contig", "")
        if isinstance(contig, list) or isinstance(contig, tuple):
            contig = "".join(contig)
        contig = self.clean(contig)
        i = 0
        answer = []
        while contig:
            if len(contig) > max_len:
                #Split lines at the commas
                pos = contig[:max_len-1].rfind(",")
                if pos == -1:
                    raise ValueError("Could not break up CONTIG")
                text, contig = contig[:pos+1], contig[pos+1:]
            else:
                text, contig = contig, ""
            answer.append(text)
        return answer

class GenBankWriter(_InsdcWriter):
    HEADER_WIDTH = 12
    QUALIFIER_INDENT = 21
    
    def _write_single_line(self, tag, text):
        "Used in the the 'header' of each GenBank record."""
        assert len(tag) < self.HEADER_WIDTH
        if len(text) > self.MAX_WIDTH - self.HEADER_WIDTH:
            import warnings
            warnings.warn("Annotation %r too long for %s line" % (text, tag))
        self.handle.write("%s%s\n" % (tag.ljust(self.HEADER_WIDTH),
                                      text.replace("\n", " ")))

    def _write_multi_line(self, tag, text):
        "Used in the the 'header' of each GenBank record."""
        #TODO - Do the line spliting while preserving white space?
        max_len = self.MAX_WIDTH - self.HEADER_WIDTH
        lines = self._split_multi_line(text, max_len)
        self._write_single_line(tag, lines[0])
        for line in lines[1:]:
            self._write_single_line("", line)

    def _write_multi_entries(self, tag, text_list):
        #used for DBLINK and any similar later line types.
        #If the list of strings is empty, nothing is written.
        for i, text in enumerate(text_list):
            if i == 0:
                self._write_single_line(tag, text)
            else:
                self._write_single_line("", text)

    def _get_date(self, record) :
        default = "01-JAN-1980"
        try :
            date = record.annotations["date"]
        except KeyError :
            return default
        #Cope with a list of one string:
        if isinstance(date, list) and len(date)==1 :
            date = date[0]
        #TODO - allow a Python date object
        if not isinstance(date, basestring) or len(date) != 11 \
        or date[2] != "-" or date[6] != "-" \
        or not date[:2].isdigit() or not date[7:].isdigit() \
        or int(date[:2]) > 31 \
        or date[3:6] not in ["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                             "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"] :
            #TODO - Check is a valid date (e.g. not 31 Feb)
            return default
        return date

    def _get_data_division(self, record):
        try:
            division = record.annotations["data_file_division"]
        except KeyError:
            division = "UNK"
        if division in ["PRI", "ROD", "MAM", "VRT", "INV", "PLN", "BCT",
                        "VRL", "PHG", "SYN", "UNA", "EST", "PAT", "STS",
                        "GSS", "HTG", "HTC", "ENV", "CON"]:
            #Good, already GenBank style
            #    PRI - primate sequences
            #    ROD - rodent sequences
            #    MAM - other mammalian sequences
            #    VRT - other vertebrate sequences
            #    INV - invertebrate sequences
            #    PLN - plant, fungal, and algal sequences
            #    BCT - bacterial sequences [plus archea]
            #    VRL - viral sequences
            #    PHG - bacteriophage sequences
            #    SYN - synthetic sequences
            #    UNA - unannotated sequences
            #    EST - EST sequences (expressed sequence tags) 
            #    PAT - patent sequences
            #    STS - STS sequences (sequence tagged sites) 
            #    GSS - GSS sequences (genome survey sequences) 
            #    HTG - HTGS sequences (high throughput genomic sequences) 
            #    HTC - HTC sequences (high throughput cDNA sequences) 
            #    ENV - Environmental sampling sequences
            #    CON - Constructed sequences
            #
            #(plus UNK for unknown)
            pass
        else:
            #See if this is in EMBL style:
            #    Division                 Code
            #    -----------------        ----
            #    Bacteriophage            PHG - common
            #    Environmental Sample     ENV - common
            #    Fungal                   FUN - map to PLN (plants + fungal)
            #    Human                    HUM - map to PRI (primates)
            #    Invertebrate             INV - common
            #    Other Mammal             MAM - common
            #    Other Vertebrate         VRT - common
            #    Mus musculus             MUS - map to ROD (rodent)
            #    Plant                    PLN - common
            #    Prokaryote               PRO - map to BCT (poor name)
            #    Other Rodent             ROD - common
            #    Synthetic                SYN - common
            #    Transgenic               TGN - ??? map to SYN ???
            #    Unclassified             UNC - map to UNK
            #    Viral                    VRL - common
            #
            #(plus XXX for submiting which we can map to UNK)
            embl_to_gbk = {"FUN":"PLN",
                           "HUM":"PRI",
                           "MUS":"ROD",
                           "PRO":"BCT",
                           "UNC":"UNK",
                           "XXX":"UNK",
                           }
            try:
                division = embl_to_gbk[division]
            except KeyError:
                division = "UNK"
        assert len(division)==3
        return division

    def _write_the_first_line(self, record):
        """Write the LOCUS line."""
        
        locus = record.name
        if not locus or locus == "<unknown name>":
            locus = record.id
        if not locus or locus == "<unknown id>":
            locus = self._get_annotation_str(record, "accession", just_first=True)
        if len(locus) > 16:
            raise ValueError("Locus identifier %r is too long" % str(locus))

        if len(record) > 99999999999:
            #Currently GenBank only officially support up to 350000, but
            #the length field can take eleven digits
            raise ValueError("Sequence too long!")

        #Get the base alphabet (underneath any Gapped or StopCodon encoding)
        a = Alphabet._get_base_alphabet(record.seq.alphabet)
        if not isinstance(a, Alphabet.Alphabet):
            raise TypeError("Invalid alphabet")
        elif isinstance(a, Alphabet.ProteinAlphabet):
            units = "aa"
        elif isinstance(a, Alphabet.NucleotideAlphabet):
            units = "bp"
        else:
            #Must be something like NucleotideAlphabet or
            #just the generic Alphabet (default for fasta files)
            raise ValueError("Need a Nucleotide or Protein alphabet")

        #Get the molecule type
        #TODO - record this explicitly in the parser?
        if isinstance(a, Alphabet.ProteinAlphabet):
            mol_type = ""
        elif isinstance(a, Alphabet.DNAAlphabet):
            mol_type = "DNA"
        elif isinstance(a, Alphabet.RNAAlphabet):
            mol_type = "RNA"
        else:
            #Must be something like NucleotideAlphabet or
            #just the generic Alphabet (default for fasta files)
            raise ValueError("Need a DNA, RNA or Protein alphabet")
        
        division = self._get_data_division(record)
        
        assert len(units) == 2
        assert len(division) == 3
        #TODO - date
        #TODO - mol_type
        line = "LOCUS       %s %s %s    %s           %s %s\n" \
                     % (locus.ljust(16),
                        str(len(record)).rjust(11),
                        units,
                        mol_type.ljust(6),
                        division,
                        self._get_date(record))
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
        assert line[55:63].strip() in ['', 'linear', 'circular'], \
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

    def _write_references(self, record):
        number = 0
        for ref in record.annotations["references"]:
            if not isinstance(ref, SeqFeature.Reference):
                continue
            number += 1
            data = str(number)
            #TODO - support more complex record reference locations?
            if ref.location and len(ref.location)==1:
                a = Alphabet._get_base_alphabet(record.seq.alphabet)
                if isinstance(a, Alphabet.ProteinAlphabet):
                    units = "residues"
                else:
                    units = "bases"
                data += "  (%s %i to %i)" % (units,
                                             ref.location[0].nofuzzy_start+1,
                                             ref.location[0].nofuzzy_end)
            self._write_single_line("REFERENCE", data)
            if ref.authors:
                #We store the AUTHORS data as a single string
                self._write_multi_line("  AUTHORS", ref.authors)
            if ref.consrtm:
                #We store the consortium as a single string
                self._write_multi_line("  CONSRTM", ref.consrtm)
            if ref.title:
                #We store the title as a single string
                self._write_multi_line("  TITLE", ref.title)
            if ref.journal:
                #We store this as a single string - holds the journal name,
                #volume, year, and page numbers of the citation
                self._write_multi_line("  JOURNAL", ref.journal)
            if ref.medline_id:
                #This line type is obsolete and was removed from the GenBank
                #flatfile format in April 2005. Should we write it?
                #Note this has a two space indent:
                self._write_multi_line("  MEDLINE", ref.medline_id)
            if ref.pubmed_id:
                #Note this has a THREE space indent:
                self._write_multi_line("   PUBMED", ref.pubmed_id)
            if ref.comment:
                self._write_multi_line("  REMARK", ref.comment)
            

    def _write_comment(self, record):
        #This is a bit complicated due to the range of possible
        #ways people might have done their annotation...
        #Currently the parser uses a single string with newlines.
        #A list of lines is also reasonable.
        #A single (long) string is perhaps the most natural of all.
        #This means we may need to deal with line wrapping.
        comment = record.annotations["comment"]
        if isinstance(comment, basestring):
            lines = comment.split("\n")
        elif isinstance(comment, list) or isinstance(comment, tuple):
            lines = comment
        else:
            raise ValueError("Could not understand comment annotation")
        self._write_multi_line("COMMENT", lines[0])
        for line in lines[1:]:
            self._write_multi_line("", line)

    def _write_contig(self, record):
        max_len = self.MAX_WIDTH - self.HEADER_WIDTH
        lines = self._split_contig(record, max_len)
        self._write_single_line("CONTIG", lines[0])
        for text in lines[1:] :
            self._write_single_line("", text)

    def _write_sequence(self, record):
        #Loosely based on code from Howard Salis
        #TODO - Force lower case?
        LETTERS_PER_LINE = 60
        SEQUENCE_INDENT = 9

        if isinstance(record.seq, UnknownSeq):
            #We have already recorded the length, and there is no need
            #to record a long sequence of NNNNNNN...NNN or whatever.
            if "contig" in record.annotations:
                self._write_contig(record)
            else:
                self.handle.write("ORIGIN\n")
            return

        #Catches sequence being None:
        data = self._get_seq_string(record).lower()
        seq_len = len(data)
        self.handle.write("ORIGIN\n")
        for line_number in range(0, seq_len, LETTERS_PER_LINE):
            self.handle.write(str(line_number+1).rjust(SEQUENCE_INDENT))
            for words in range(line_number,
                               min(line_number+LETTERS_PER_LINE, seq_len), 10):
                self.handle.write(" %s" % data[words:words+10])
            self.handle.write("\n")
        
    def write_record(self, record):
        """Write a single record to the output file."""
        handle = self.handle
        self._write_the_first_line(record)

        accession = self._get_annotation_str(record, "accession",
                                             record.id.split(".", 1)[0],
                                             just_first=True)
        acc_with_version = accession
        if record.id.startswith(accession+"."):
            try:
                acc_with_version = "%s.%i" \
                                   % (accession,
                                      int(record.id.split(".", 1)[1]))
            except ValueError:
                pass
        gi = self._get_annotation_str(record, "gi", just_first=True)

        descr = record.description
        if descr == "<unknown description>" : descr = "."
        self._write_multi_line("DEFINITION", descr)
        
        self._write_single_line("ACCESSION", accession)
        if gi != ".":
            self._write_single_line("VERSION", "%s  GI:%s" \
                                    % (acc_with_version, gi))
        else:
            self._write_single_line("VERSION", "%s" % (acc_with_version))

        #The NCBI only expect two types of link so far,
        #e.g. "Project:28471" and "Trace Assembly Archive:123456"
        #TODO - Filter the dbxrefs list to just these?
        self._write_multi_entries("DBLINK", record.dbxrefs)

        try:
            #List of strings
            #Keywords should be given separated with semi colons,
            keywords = "; ".join(record.annotations["keywords"])
            #with a trailing period:
            if not keywords.endswith(".") :
                keywords += "."
        except KeyError:
            #If no keywords, there should be just a period:
            keywords = "."
        self._write_multi_line("KEYWORDS", keywords)

        if "segment" in record.annotations:
            #Deal with SEGMENT line found only in segmented records,
            #e.g. AH000819
            segment = record.annotations["segment"]
            if isinstance(segment, list):
                assert len(segment)==1, segment
                segment = segment[0]
            self._write_single_line("SEGMENT", segment)

        self._write_multi_line("SOURCE", \
                                self._get_annotation_str(record, "source"))
        #The ORGANISM line MUST be a single line, as any continuation is the taxonomy
        org = self._get_annotation_str(record, "organism")
        if len(org) > self.MAX_WIDTH - self.HEADER_WIDTH:
            org = org[:self.MAX_WIDTH - self.HEADER_WIDTH-4]+"..."
        self._write_single_line("  ORGANISM", org)
        try:
            #List of strings
            #Taxonomy should be given separated with semi colons,
            taxonomy = "; ".join(record.annotations["taxonomy"])
            #with a trailing period:
            if not taxonomy.endswith(".") :
                taxonomy += "."
        except KeyError:
            taxonomy = "."
        self._write_multi_line("", taxonomy)

        if "references" in record.annotations:
            self._write_references(record)

        if "comment" in record.annotations:
            self._write_comment(record)

        handle.write("FEATURES             Location/Qualifiers\n")
        rec_length = len(record)
        for feature in record.features:
            self._write_feature(feature, rec_length) 
        self._write_sequence(record)
        handle.write("//\n")

class EmblWriter(_InsdcWriter):
    HEADER_WIDTH = 5
    QUALIFIER_INDENT = 21
    QUALIFIER_INDENT_STR = "FT" + " "*(QUALIFIER_INDENT-2)
    QUALIFIER_INDENT_TMP = "FT   %s                " # 21 if %s is empty
    FEATURE_HEADER = "FH   Key             Location/Qualifiers\n"
    
    def _write_contig(self, record):
        max_len = self.MAX_WIDTH - self.HEADER_WIDTH
        lines = self._split_contig(record, max_len)
        for text in lines:
            self._write_single_line("CO", text)

    def _write_sequence(self, record):
        LETTERS_PER_BLOCK = 10
        BLOCKS_PER_LINE = 6
        LETTERS_PER_LINE = LETTERS_PER_BLOCK * BLOCKS_PER_LINE
        POSITION_PADDING = 10
        handle = self.handle #save looking up this multiple times
        
        if isinstance(record.seq, UnknownSeq):
            #We have already recorded the length, and there is no need
            #to record a long sequence of NNNNNNN...NNN or whatever.
            if "contig" in record.annotations:
                self._write_contig(record)
            else:
                #TODO - Can the sequence just be left out as in GenBank files?
                handle.write("SQ   \n")
            return

        #Catches sequence being None
        data = self._get_seq_string(record).lower()
        seq_len = len(data)

        #Get the base alphabet (underneath any Gapped or StopCodon encoding)
        a = Alphabet._get_base_alphabet(record.seq.alphabet)
        if isinstance(a, Alphabet.DNAAlphabet):
            #TODO - What if we have RNA?
            a_count = data.count('A') + data.count('a')
            c_count = data.count('C') + data.count('c')
            g_count = data.count('G') + data.count('g')
            t_count = data.count('T') + data.count('t')
            other = seq_len - (a_count + c_count + g_count + t_count)
            handle.write("SQ   Sequence %i BP; %i A; %i C; %i G; %i T; %i other;\n" \
                         % (seq_len, a_count, c_count, g_count, t_count, other))
        else:
            handle.write("SQ   \n")
        
        for line_number in range(0, seq_len // LETTERS_PER_LINE):
            handle.write("    ") #Just four, not five
            for block in range(BLOCKS_PER_LINE) :
                index = LETTERS_PER_LINE*line_number + LETTERS_PER_BLOCK*block
                handle.write((" %s" % data[index:index+LETTERS_PER_BLOCK]))
            handle.write(str((line_number+1)
                             *LETTERS_PER_LINE).rjust(POSITION_PADDING))
            handle.write("\n")
        if seq_len % LETTERS_PER_LINE:
            #Final (partial) line
            line_number = (seq_len // LETTERS_PER_LINE)
            handle.write("    ") #Just four, not five
            for block in range(BLOCKS_PER_LINE) :
                index = LETTERS_PER_LINE*line_number + LETTERS_PER_BLOCK*block
                handle.write((" %s" % data[index:index+LETTERS_PER_BLOCK]).ljust(11))
            handle.write(str(seq_len).rjust(POSITION_PADDING))
            handle.write("\n")

    def _write_single_line(self, tag, text):
        assert len(tag)==2
        line = tag+"   "+text
        if len(text) > self.MAX_WIDTH:
            import warnings
            warnings.warn("Line %r too long" % line)
        self.handle.write(line+"\n")

    def _write_multi_line(self, tag, text):
        max_len = self.MAX_WIDTH - self.HEADER_WIDTH
        lines = self._split_multi_line(text, max_len)
        for line in lines :
            self._write_single_line(tag, line)
        
    def _write_the_first_lines(self, record):
        """Write the ID and AC lines."""
        if "." in record.id and record.id.rsplit(".", 1)[1].isdigit():
            version = "SV " + record.id.rsplit(".", 1)[1]
            accession = self._get_annotation_str(record, "accession",
                                                 record.id.rsplit(".", 1)[0],
                                                 just_first=True)
        else :
            version = ""
            accession = self._get_annotation_str(record, "accession",
                                                 record.id,
                                                 just_first=True)
        
        if ";" in accession :
            raise ValueError("Cannot have semi-colon in EMBL accession, %s" \
                             % repr(str(accession)))
        if " " in accession :
            #This is out of practicallity... might it be allowed?
            raise ValueError("Cannot have spaces in EMBL accession, %s" \
                             % repr(str(accession)))

        #Get the molecule type
        #TODO - record this explicitly in the parser?
        #Get the base alphabet (underneath any Gapped or StopCodon encoding)
        a = Alphabet._get_base_alphabet(record.seq.alphabet)
        if not isinstance(a, Alphabet.Alphabet):
            raise TypeError("Invalid alphabet")
        elif isinstance(a, Alphabet.DNAAlphabet):
            mol_type = "DNA"
            units = "BP"
        elif isinstance(a, Alphabet.RNAAlphabet):
            mol_type = "RNA"
            units = "BP"
        elif isinstance(a, Alphabet.ProteinAlphabet):
            mol_type = "PROTEIN"
            units = "AA"
        else:
            #Must be something like NucleotideAlphabet
            raise ValueError("Need a DNA, RNA or Protein alphabet")

        #Get the taxonomy division
        division = self._get_data_division(record)

        #TODO - Full ID line
        handle = self.handle
        #ID   <1>; SV <2>; <3>; <4>; <5>; <6>; <7> BP.
        #1. Primary accession number
        #2. Sequence version number
        #3. Topology: 'circular' or 'linear'
        #4. Molecule type
        #5. Data class
        #6. Taxonomic division
        #7. Sequence length
        self._write_single_line("ID", "%s; %s; ; %s; ; %s; %i %s." \
                                % (accession, version, mol_type,
                                   division, len(record), units))
        handle.write("XX\n")
        self._write_single_line("AC", accession+";")
        handle.write("XX\n")

    def _get_data_division(self, record):
        try:
            division = record.annotations["data_file_division"]
        except KeyError:
            division = "UNC"
        if division in ["PHG", "ENV", "FUN", "HUM", "INV", "MAM", "VRT",
                        "MUS", "PLN", "PRO", "ROD", "SYN", "TGN", "UNC",
                        "VRL", "XXX"]:
            #Good, already EMBL style
            #    Division                 Code
            #    -----------------        ----
            #    Bacteriophage            PHG
            #    Environmental Sample     ENV
            #    Fungal                   FUN
            #    Human                    HUM
            #    Invertebrate             INV
            #    Other Mammal             MAM
            #    Other Vertebrate         VRT
            #    Mus musculus             MUS
            #    Plant                    PLN
            #    Prokaryote               PRO
            #    Other Rodent             ROD
            #    Synthetic                SYN
            #    Transgenic               TGN
            #    Unclassified             UNC (i.e. unknown)
            #    Viral                    VRL
            #
            #(plus XXX used for submiting data to EMBL)
            pass
        else:
            #See if this is in GenBank style & can be converted.
            #Generally a problem as the GenBank groups are wider
            #than those of EMBL. Note that GenBank use "BCT" for
            #both bacteria and acherea thus this maps to EMBL's
            #"PRO" nicely.
            gbk_to_embl = {"BCT":"PRO",
                           "UNK":"UNC",
                           }
            try:
                division = gbk_to_embl[division]
            except KeyError:
                division = "UNC"
        assert len(division)==3
        return division

    def _write_references(self, record):
        #The order should be RN, RC, RP, RX, RG, RA, RT, RL
        number = 0
        for ref in record.annotations["references"]:
            if not isinstance(ref, SeqFeature.Reference):
                continue
            number += 1
            self._write_single_line("RN", "[%i]" % number)
            #TODO - support for RC line (needed in parser too)
            #TODO - support more complex record reference locations?
            if ref.location and len(ref.location)==1:
                self._write_single_line("RP", "%i-%i" % (ref.location[0].nofuzzy_start+1,
                                                         ref.location[0].nofuzzy_end))
            #TODO - record any DOI or AGRICOLA identifier in the reference object?
            if ref.pubmed_id:
                self._write_single_line("RX", "PUBMED; %s." % ref.pubmed_id)
            if ref.consrtm:
                self._write_single_line("RG", "%s" % ref.consrtm)
            if ref.authors:
                #We store the AUTHORS data as a single string
                self._write_multi_line("RA", ref.authors+";")
            if ref.title:
                #We store the title as a single string
                self._write_multi_line("RT", '"%s";' % ref.title)
            if ref.journal:
                #We store this as a single string - holds the journal name,
                #volume, year, and page numbers of the citation
                self._write_multi_line("RL", ref.journal)
            self.handle.write("XX\n")

    def _write_comment(self, record):
        #This is a bit complicated due to the range of possible
        #ways people might have done their annotation...
        #Currently the parser uses a single string with newlines.
        #A list of lines is also reasonable.
        #A single (long) string is perhaps the most natural of all.
        #This means we may need to deal with line wrapping.
        comment = record.annotations["comment"]
        if isinstance(comment, basestring):
            lines = comment.split("\n")
        elif isinstance(comment, list) or isinstance(comment, tuple):
            lines = comment
        else:
            raise ValueError("Could not understand comment annotation")
        #TODO - Merge this with the GenBank comment code?
        if not lines : return
        for line in lines:
            self._write_multi_line("CC", line)
        self.handle.write("XX\n")

    def write_record(self, record):
        """Write a single record to the output file."""

        handle = self.handle
        self._write_the_first_lines(record)

        #PR line (0 or 1 lines only), project identifier
        for xref in record.dbxrefs:
            if xref.startswith("Project:"):
                self._write_single_line("PR", xref+";")
                handle.write("XX\n")
                break

        #TODO - DT lines (date)

        descr = record.description
        if descr == "<unknown description>" : descr = "."
        self._write_multi_line("DE", descr)
        handle.write("XX\n")

        #Should this be "source" or "organism"?
        self._write_multi_line("OS", self._get_annotation_str(record, "organism"))
        try:
            #List of strings
            taxonomy = "; ".join(record.annotations["taxonomy"]) + "."
        except KeyError:
            taxonomy = "."
        self._write_multi_line("OC", taxonomy)
        handle.write("XX\n")

        if "references" in record.annotations:
            self._write_references(record)

        if "comment" in record.annotations:
            self._write_comment(record)

        handle.write(self.FEATURE_HEADER)
        rec_length = len(record)
        for feature in record.features:
            self._write_feature(feature, rec_length)

        self._write_sequence(record)
        handle.write("//\n")

class ImgtWriter(EmblWriter):
    HEADER_WIDTH = 5
    QUALIFIER_INDENT = 25 # Not 21 as in EMBL
    QUALIFIER_INDENT_STR = "FT" + " "*(QUALIFIER_INDENT-2)
    QUALIFIER_INDENT_TMP = "FT   %s                    " # 25 if %s is empty
    FEATURE_HEADER = "FH   Key                 Location/Qualifiers\n"

if __name__ == "__main__":
    print "Quick self test"
    import os
    from StringIO import StringIO

    def compare_record(old, new):
        if old.id != new.id and old.name != new.name:
            raise ValueError("'%s' or '%s' vs '%s' or '%s' records" \
                             % (old.id, old.name, new.id, new.name))
        if len(old.seq) != len(new.seq):
            raise ValueError("%i vs %i" % (len(old.seq), len(new.seq)))
        if str(old.seq).upper() != str(new.seq).upper():
            if len(old.seq) < 200:
                raise ValueError("'%s' vs '%s'" % (old.seq, new.seq))
            else:
                raise ValueError("'%s...' vs '%s...'" % (old.seq[:100], new.seq[:100]))
        if old.features and new.features:
            return compare_features(old.features, new.features)
        #Just insist on at least one word in common:
        if (old.description or new.description) \
        and not set(old.description.split()).intersection(new.description.split()):
            raise ValueError("%s versus %s" \
                             % (repr(old.description), repr(new.description)))
        #TODO - check annotation
        if "contig" in old.annotations:
            assert old.annotations["contig"] == \
                   new.annotations["contig"]
        return True

    def compare_records(old_list, new_list):
        """Check two lists of SeqRecords agree, raises a ValueError if mismatch."""
        if len(old_list) != len(new_list):
            raise ValueError("%i vs %i records" % (len(old_list), len(new_list)))
        for old, new in zip(old_list, new_list):
            if not compare_record(old, new):
                return False
        return True

    def compare_feature(old, new, ignore_sub_features=False):
        """Check two SeqFeatures agree."""
        if old.type != new.type:
            raise ValueError("Type %s versus %s" % (old.type, new.type))
        if old.location.nofuzzy_start != new.location.nofuzzy_start \
        or old.location.nofuzzy_end != new.location.nofuzzy_end:
            raise ValueError("%s versus %s:\n%s\nvs:\n%s" \
                             % (old.location, new.location, str(old), str(new)))
        if old.strand != new.strand:
            raise ValueError("Different strand:\n%s\nvs:\n%s" % (str(old), str(new)))
        if old.location.start != new.location.start:
            raise ValueError("Start %s versus %s:\n%s\nvs:\n%s" \
                             % (old.location.start, new.location.start, str(old), str(new)))
        if old.location.end != new.location.end:
            raise ValueError("End %s versus %s:\n%s\nvs:\n%s" \
                             % (old.location.end, new.location.end, str(old), str(new)))
        if not ignore_sub_features:
            if len(old.sub_features) != len(new.sub_features):
                raise ValueError("Different sub features")
            for a, b in zip(old.sub_features, new.sub_features):
                if not compare_feature(a, b):
                    return False
        #This only checks key shared qualifiers
        #Would a white list be easier?
        #for key in ["name", "gene", "translation", "codon_table", "codon_start", "locus_tag"]:
        for key in set(old.qualifiers).intersection(new.qualifiers):
            if key in ["db_xref", "protein_id", "product", "note"]:
                #EMBL and GenBank files are use different references/notes/etc
                continue
            if old.qualifiers[key] != new.qualifiers[key]:
                raise ValueError("Qualifier mis-match for %s:\n%s\n%s" \
                                 % (key, old.qualifiers[key], new.qualifiers[key]))
        return True

    def compare_features(old_list, new_list, ignore_sub_features=False):
        """Check two lists of SeqFeatures agree, raises a ValueError if mismatch."""
        if len(old_list) != len(new_list):
            raise ValueError("%i vs %i features" % (len(old_list), len(new_list)))
        for old, new in zip(old_list, new_list):
            #This assumes they are in the same order
            if not compare_feature(old, new, ignore_sub_features):
                return False
        return True

    def check_genbank_writer(records):
        handle = StringIO()
        GenBankWriter(handle).write_file(records)
        handle.seek(0)

        records2 = list(GenBankIterator(handle))
        assert compare_records(records, records2)

    def check_embl_writer(records):
        handle = StringIO()
        try:
            EmblWriter(handle).write_file(records)
        except ValueError, err:
            print err
            return
        handle.seek(0)

        records2 = list(EmblIterator(handle))
        assert compare_records(records, records2)

    for filename in os.listdir("../../Tests/GenBank"):
        if not filename.endswith(".gbk") and not filename.endswith(".gb"):
            continue
        print filename
        
        handle = open("../../Tests/GenBank/%s" % filename)
        records = list(GenBankIterator(handle))
        handle.close()

        check_genbank_writer(records)
        check_embl_writer(records)

    for filename in os.listdir("../../Tests/EMBL"):
        if not filename.endswith(".embl"):
            continue
        print filename
        
        handle = open("../../Tests/EMBL/%s" % filename)
        records = list(EmblIterator(handle))
        handle.close()

        check_genbank_writer(records)
        check_embl_writer(records)

    from Bio import SeqIO
    for filename in os.listdir("../../Tests/SwissProt"):
        if not filename.startswith("sp"):
            continue
        print filename
        
        handle = open("../../Tests/SwissProt/%s" % filename)
        records = list(SeqIO.parse(handle, "swiss"))
        handle.close()

        check_genbank_writer(records)

