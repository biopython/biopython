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
from Bio import SeqFeature

# NOTE
# ====
# The "brains" for parsing GenBank and EMBL files (and any
# other flat file variants from the INSDC in future) is in
# Bio.GenBank.Scanner (plus the _FeatureConsumer in Bio.GenBank)
# However, all the writing code is in this file.


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

def _insdc_feature_position_string(pos, offset=0):
    """Build a GenBank/EMBL position string (PRIVATE).

    Use offset=1 to add one to convert a start position from python counting.
    """
    if isinstance(pos, SeqFeature.ExactPosition) :
        return "%i" % (pos.position+offset)
    elif isinstance(pos, SeqFeature.WithinPosition) :
        return "(%i.%i)" % (pos.position + offset,
                            pos.position + pos.extension + offset)
    elif isinstance(pos, SeqFeature.BetweenPosition) :
        return "(%i^%i)" % (pos.position + offset,
                            pos.position + pos.extension + offset)
    elif isinstance(pos, SeqFeature.BeforePosition) :
        return "<%i" % (pos.position + offset)
    elif isinstance(pos, SeqFeature.AfterPosition) :
        return ">%i" % (pos.position + offset)
    elif isinstance(pos, SeqFeature.OneOfPosition):
        return "one-of(%s)" \
               % ",".join([_insdc_feature_position_string(p,offset) \
                           for p in pos.position_choices])
    elif isinstance(pos, SeqFeature.AbstractPosition) :
        raise NotImplementedError("Please report this as a bug in Biopython.")
    else :
        raise ValueError("Expected a SeqFeature position object.")


def _insdc_location_string_ignoring_strand_and_subfeatures(feature) :
    if feature.ref :
        ref = "%s:" % feature.ref
    else :
        ref = ""
    assert not feature.ref_db
    if feature.location.start==feature.location.end \
    and isinstance(feature.location.end, SeqFeature.ExactPosition):
        #Special case, 12^13 gets mapped to location 12:12
        #(a zero length slice, meaning the point between two letters)
        return "%s%i^%i" % (ref, feature.location.end.position,
                            feature.location.end.position+1)
    else :
        #Typical case, e.g. 12..15 gets mapped to 11:15
        return ref \
               + _insdc_feature_position_string(feature.location.start, +1) \
               + ".." + \
               _insdc_feature_position_string(feature.location.end)

def _insdc_feature_location_string(feature):
    """Build a GenBank/EMBL location string from a SeqFeature (PRIVATE)."""
    # Have a choice of how to show joins on the reverse complement strand,
    # complement(join(1,10),(20,100)) vs join(complement(20,100),complement(1,10))
    # Notice that the order of the entries gets flipped!
    #
    # GenBank and EMBL would both use now complement(join(1,10),(20,100))
    # which is shorter at least.
    #
    # In the above situations, we expect the parent feature and the two children
    # to all be marked as strand==-1, and in the order 0:10 then 19:100.
    #
    # Also need to consider dual-strand examples like these from the Arabidopsis
    # thaliana chloroplast NC_000932: join(complement(69611..69724),139856..140650)
    # gene ArthCp047, GeneID:844801 or its CDS which is even better due to a splice:
    # join(complement(69611..69724),139856..140087,140625..140650)
    # protein NP_051038.1 GI:7525057
    #

    if not feature.sub_features :
        #Non-recursive.
        #assert feature.location_operator == "", \
        #       "%s has no subfeatures but location_operator %s" \
        #       % (repr(feature), feature.location_operator)
        location = _insdc_location_string_ignoring_strand_and_subfeatures(feature)
        if feature.strand == -1 :
            location = "complement(%s)" % location
        return location
    # As noted above, treat reverse complement strand features carefully:
    if feature.strand == -1 :
        for f in feature.sub_features :
            assert f.strand == -1
        return "complement(%s(%s))" \
               % (feature.location_operator,
                  ",".join(_insdc_location_string_ignoring_strand_and_subfeatures(f) \
                           for f in feature.sub_features))
    #if feature.strand == +1 :
    #    for f in feature.sub_features :
    #        assert f.strand == +1
    #This covers typical forward strand features, and also an evil mixed strand:
    assert feature.location_operator != ""
    return  "%s(%s)" % (feature.location_operator,
                        ",".join([_insdc_feature_location_string(f) \
                                  for f in feature.sub_features]))


class GenBankWriter(SequentialSequenceWriter) :
    HEADER_WIDTH = 12
    MAX_WIDTH = 80
    QUALIFIER_INDENT = 21
    
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

    def _write_multi_entries(self, tag, text_list) :
        #used for DBLINK and any similar later line types.
        #If the list of strings is empty, nothing is written.
        for i, text in enumerate(text_list) :
            if i==0 :
                self._write_single_line(tag, text)
            else :
                self._write_single_line("", text)

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
            units = "aa"
        elif isinstance(a, Alphabet.NucleotideAlphabet) :
            units = "bp"
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
                            "GSS","HTG","HTC","ENV","CON"] :
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

        #The NCBI only expect two types of link so far,
        #e.g. "Project:28471" and "Trace Assembly Archive:123456"
        #TODO - Filter the dbxrefs list to just these?
        self._write_multi_entries("DBLINK", record.dbxrefs)

        try :
            #List of strings
            keywords = "; ".join(record.annotations["keywords"])
        except KeyError :
            keywords = "."
        self._write_multi_line("KEYWORDS", keywords)

        if "segment" in record.annotations :
            #Deal with SEGMENT line found only in segmented records,
            #e.g. AH000819
            segment = record.annotations["segment"]
            if isinstance(segment, list) :
                assert len(segment)==1, segment
                segment = segment[0]
            self._write_single_line("SEGMENT", segment)

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

    def _write_feature_qualifier(self, key, value=None, quote=None) :
        if not value :
            self.handle.write("%s/%s\n" % (" "*self.QUALIFIER_INDENT, key))
            return
        #Quick hack with no line wrapping, may be useful for testing:
        #self.handle.write('%s/%s="%s"\n' % (" "*self.QUALIFIER_INDENT, key, value))
        if quote is None :
            #Try to mimic unwritten rules about when quotes can be left out:
            if isinstance(value, int) or isinstance(value, long) :
                quote = False
            else :
                quote = True
        if quote :
            line = '%s/%s="%s"' % (" "*self.QUALIFIER_INDENT, key, value)
        else :
            line = '%s/%s=%s' % (" "*self.QUALIFIER_INDENT, key, value)
        if len(line) < self.MAX_WIDTH :
            self.handle.write(line+"\n")
            return
        while line.lstrip() :
            if len(line) < self.MAX_WIDTH :
                self.handle.write(line+"\n")
                return
            #Insert line break...
            for index in range(min(len(line)-1,self.MAX_WIDTH),self.QUALIFIER_INDENT+1,-1) :
                if line[index]==" " : break
            if line[index] != " " :
                #No nice place to break...
                index = self.MAX_WIDTH
            self.handle.write(line[:index] + "\n")
            line = " "*self.QUALIFIER_INDENT + line[index:].lstrip()

    def _wrap_location(self, location) :
        """Split a feature location into lines (break at commas)."""
        #TODO - Rewrite this not to recurse!
        length = self.MAX_WIDTH - self.QUALIFIER_INDENT
        if len(location) <= length :
            return location
        index = location[:length].rfind(",")
        if index == -1 :
            #No good place to split (!)
            import warnings
            warnings.warn("Couldn't split location:\n%s" % location)
            return location
        return location[:index+1] + "\n" + \
               " "*self.QUALIFIER_INDENT + self._wrap_location(location[index+1:])

    def _write_feature(self, feature):
        """Write a single SeqFeature object to features table."""
        assert feature.type, feature
        #TODO - Line wrapping for long locations!
        location = _insdc_feature_location_string(feature)
        line = ("     %s                " % feature.type)[:self.QUALIFIER_INDENT] \
               + self._wrap_location(location) + "\n"
        self.handle.write(line)
        #Now the qualifiers...
        for key, values in feature.qualifiers.iteritems() :
            if isinstance(values, list) or isinstance(values, tuple) :
                for value in values :
                    self._write_feature_qualifier(key, value)
            elif values :
                #String, int, etc
                self._write_feature_qualifier(key, values)
            else :
                #e.g. a /psuedo entry
                self._write_feature_qualifier(key)


if __name__ == "__main__" :
    print "Quick self test"
    import os
    from StringIO import StringIO

    def compare_record(old, new) :
        if old.id != new.id and old.name != new.name :
            raise ValueError("'%s' or '%s' vs '%s' or '%s' records" \
                             % (old.id, old.name, new.id, new.name))
        if len(old.seq) != len(new.seq) :
            raise ValueError("%i vs %i" % (len(old.seq), len(new.seq)))
        if str(old.seq).upper() != str(new.seq).upper() :
            if len(old.seq) < 200 :
                raise ValueError("'%s' vs '%s'" % (old.seq, new.seq))
            else :
                raise ValueError("'%s...' vs '%s...'" % (old.seq[:100], new.seq[:100]))
        if old.features and new.features :
            return compare_features(old.features, new.features)
        #Just insist on at least one word in common:
        if (old.description or new.description) \
        and not set(old.description.split()).intersection(new.description.split()):
            raise ValueError("%s versus %s" \
                             % (repr(old.description), repr(new.description)))
        #TODO - check annotation
        return True

    def compare_records(old_list, new_list) :
        """Check two lists of SeqRecords agree, raises a ValueError if mismatch."""
        if len(old_list) != len(new_list) :
            raise ValueError("%i vs %i records" % (len(old_list), len(new_list)))
        for old, new in zip(old_list, new_list) :
            if not compare_record(old,new) :
                return False
        return True

    def compare_feature(old, new, ignore_sub_features=False) :
        """Check two SeqFeatures agree."""
        if old.type != new.type :
            raise ValueError("Type %s versus %s" % (old.type, new.type))
        if old.location.nofuzzy_start != new.location.nofuzzy_start \
        or old.location.nofuzzy_end != new.location.nofuzzy_end :
            raise ValueError("%s versus %s:\n%s\nvs:\n%s" \
                             % (old.location, new.location, str(old), str(new)))
        if old.strand != new.strand :
            raise ValueError("Different strand:\n%s\nvs:\n%s" % (str(old), str(new)))
        if old.location.start != new.location.start :
            raise ValueError("Start %s versus %s:\n%s\nvs:\n%s" \
                             % (old.location.start, new.location.start, str(old), str(new)))
        if old.location.end != new.location.end :
            raise ValueError("End %s versus %s:\n%s\nvs:\n%s" \
                             % (old.location.end, new.location.end, str(old), str(new)))
        if not ignore_sub_features :
            if len(old.sub_features) != len(new.sub_features) :
                raise ValueError("Different sub features")
            for a,b in zip(old.sub_features, new.sub_features) :
                if not compare_feature(a,b) :
                    return False
        #This only checks key shared qualifiers
        #Would a white list be easier?
        #for key in ["name","gene","translation","codon_table","codon_start","locus_tag"] :
        for key in set(old.qualifiers.keys()).intersection(new.qualifiers.keys()):
            if key in ["db_xref","protein_id","product","note"] :
                #EMBL and GenBank files are use different references/notes/etc
                continue
            if old.qualifiers[key] != new.qualifiers[key] :
                raise ValueError("Qualifier mis-match for %s:\n%s\n%s" \
                                 % (key, old.qualifiers[key], new.qualifiers[key]))
        return True

    def compare_features(old_list, new_list, ignore_sub_features=False) :
        """Check two lists of SeqFeatures agree, raises a ValueError if mismatch."""
        if len(old_list) != len(new_list) :
            raise ValueError("%i vs %i features" % (len(old_list), len(new_list)))
        for old, new in zip(old_list, new_list) :
            #This assumes they are in the same order
            if not compare_feature(old,new,ignore_sub_features) :
                return False
        return True

    def check_genbank_writer(records) :
        handle = StringIO()
        GenBankWriter(handle).write_file(records)
        handle.seek(0)

        records2 = list(GenBankIterator(handle))
        assert compare_records(records, records2)

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

