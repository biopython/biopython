# Copyright 2006-2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.Align.Generic import Alignment
from Interfaces import AlignmentIterator, SequentialAlignmentWriter

class StockholmWriter(SequentialAlignmentWriter) :
    """Stockholm/PFAM alignment writer."""

    #These dictionaries should be kept in sync with those
    #defined in the StockholmIterator class.
    pfam_gr_mapping = {"secondary_structure" : "SS",
                       "surface_accessibility" : "SA",
                       "transmembrane" : "TM",
                       "posterior_probability" : "PP",
                       "ligand_binding" : "LI",
                       "active_site" : "AS",
                       "intron" : "IN"}
    #Following dictionary deliberately does not cover AC, DE or DR
    pfam_gs_mapping = {"organism" : "OS",
                       "organism_classification" : "OC",
                       "look" : "LO"}

    def write_alignment(self, alignment) :
        """Use this to write (another) single alignment to an open file.
        
        Note that sequences and their annotation are recorded
        together (rather than having a block of annotation followed
        by a block of aligned sequences).
        """
        records = alignment.get_all_seqs()
        count = len(records)
        
        self._length_of_sequences = alignment.get_alignment_length()
        self._ids_written = []
        
        #NOTE - For now, the alignment object does not hold any per column
        #or per alignment annotation - only per sequence.
        
        if count == 0 :
            raise ValueError("Must have at least one sequence")
        if self._length_of_sequences == 0 :
            raise ValueError("Non-empty sequences are required")

        self.handle.write("# STOCKHOLM 1.0\n")
        self.handle.write("#=GF SQ %i\n" % count)
        for record in records :
            self._write_record(record)
        self.handle.write("//\n")

    def _write_record(self, record):
        """Write a single SeqRecord to the file"""
        if self._length_of_sequences <> len(record.seq) :
            raise ValueError("Sequences must all be the same length")

        #For the case for stockholm to stockholm, try and use record.name
        seq_name = record.id
        if record.name is not None :
            if "accession" in record.annotations :
                if record.id == record.annotations["accession"] :
                    seq_name = record.name

        #In the Stockholm file format, spaces are not allowed in the id
        seq_name = seq_name.replace(" ","_")

        if "start" in record.annotations \
        and  "end" in record.annotations :
            suffix = "/%s-%s" % (str(record.annotations["start"]),
                                 str(record.annotations["end"]))
            if seq_name[-len(suffix):] <> suffix :
                seq_name = "%s/%s-%s" % (seq_name,
                                        str(record.annotations["start"]),
                                        str(record.annotations["end"]))

        if seq_name in self._ids_written :
            raise ValueError("Duplicate record identifier: %s" % seq_name)
        self._ids_written.append(seq_name)
        self.handle.write("%s %s\n" % (seq_name, record.seq.tostring()))

        #The recommended placement for GS lines (per sequence annotation)
        #is above the alignment (as a header block) or just below the
        #corresponding sequence.
        #
        #The recommended placement for GR lines (per sequence per column
        #annotation such as secondary structure) is just below the
        #corresponding sequence.
        #
        #We put both just below the corresponding sequence as this allows
        #us to write the file using a single pass through the records.

        #AC = Accession
        if "accession" in record.annotations :
            self.handle.write("#=GS %s AC %s\n" \
                % (seq_name, self.clean(record.annotations["accession"])))
        elif record.id :
            self.handle.write("#=GS %s AC %s\n" \
                % (seq_name, self.clean(record.id)))
        
        #DE = description
        if record.description :
            self.handle.write("#=GS %s DE %s\n" \
                % (seq_name, self.clean(record.description)))

        #DE = database links
        for xref in record.dbxrefs :
            self.handle.write("#=GS %s DR %s\n" \
                % (seq_name, self.clean(xref)))

        #GS/GR = other per sequence annotation
        for key in record.annotations :
            if key in self.pfam_gs_mapping :
                self.handle.write("#=GS %s %s %s\n" \
                                  % (seq_name,
                                     self.clean(self.pfam_gs_mapping[key]),
                                     self.clean(str(record.annotations[key]))))
            elif key in self.pfam_gr_mapping :
                if len(str(record.annotations[key]))==len(record.seq) :
                    self.handle.write("#=GR %s %s %s\n" \
                                      % (seq_name,
                                         self.clean(self.pfam_gr_mapping[key]),
                                         self.clean((record.annotations[key]))))
                else :
                    #Should we print a warning?
                    pass
            else :
                #It doesn't follow the PFAM standards, but should we record
                #this data anyway?
                pass
        
class StockholmIterator(AlignmentIterator) :
    """Loads a Stockholm file from PFAM into Alignment objects.

    The file may contain multiple concatenated alignments, which are loaded
    and returned incrementally.

    This parser will detect if the Stockholm file follows the PFAM
    conventions for sequence specific meta-data (lines starting #=GS
    and #=GR) and populates the SeqRecord fields accordingly.
    
    Any annotation which does not follow the PFAM conventions is currently
    ignored.

    If an accession is provided for an entry in the meta data, IT WILL NOT
    be used as the record.id (it will be recorded in the record's
    annotations).  This is because some files have (sub) sequences from
    different parts of the same accession (differentiated by different
    start-end positions).

    Wrap-around alignments are not supported - each sequences must be on
    a single line.  However, interlaced sequences should work.

    For more information on the file format, please see:
    http://www.bioperl.org/wiki/Stockholm_multiple_alignment_format
    http://www.cgb.ki.se/cgb/groups/sonnhammer/Stockholm.html

    For consistency with BioPerl and EMBOSS we call this the "stockholm"
    format.
    """

    #These dictionaries should be kept in sync with those
    #defined in the PfamStockholmWriter class.
    pfam_gr_mapping = {"SS" : "secondary_structure",
                       "SA" : "surface_accessibility",
                       "TM" : "transmembrane",
                       "PP" : "posterior_probability",
                       "LI" : "ligand_binding",
                       "AS" : "active_site",
                       "IN" : "intron"}
    #Following dictionary deliberately does not cover AC, DE or DR
    pfam_gs_mapping = {"OS" : "organism",
                       "OC" : "organism_classification",
                       "LO" : "look"}

    def next(self) :
        try :
            line = self._header
            del self._header
        except AttributeError :
            line = self.handle.readline()
        if not line:
            #Empty file - just give up.
            return
        if not line.strip() == '# STOCKHOLM 1.0':
            raise ValueError("Did not find STOCKHOLM header")
            #import sys
            #print >> sys.stderr, 'Warning file does not start with STOCKHOLM 1.0'

        # Note: If this file follows the PFAM conventions, there should be
        # a line containing the number of sequences, e.g. "#=GF SQ 67"
        # We do not check for this - perhaps we should, and verify that
        # if present it agrees with our parsing.

        seqs = {}
        ids = []
        gs = {}
        gr = {}
        gf = {}
        passed_end_alignment = False
        while 1:
            line = self.handle.readline()
            if not line: break #end of file
            line = line.strip() #remove trailing \n
            if line == '# STOCKHOLM 1.0':
                self._header = line
                break
            elif line == "//" :
                #The "//" line indicates the end of the alignment.
                #There may still be more meta-data
                passed_end_alignment = True
            elif line == "" :
                #blank line, ignore
                pass
            elif line[0] <> "#" :
                #Sequence
                #Format: "<seqname> <sequence>"
                assert not passed_end_alignment
                parts = [x.strip() for x in line.split(" ",1)]
                if len(parts) <> 2 :
                    #This might be someone attempting to store a zero length sequence?
                    raise ValueError("Could not split line into identifier " \
                                      + "and sequence:\n" + line)
                id, seq = parts
                if id not in ids :
                    ids.append(id)
                seqs.setdefault(id, '')
                seqs[id] += seq.replace(".","-")
            elif len(line) >= 5 :
                #Comment line or meta-data
                if line[:5] == "#=GF " :
                    #Generic per-File annotation, free text
                    #Format: #=GF <feature> <free text>
                    feature, text = line[5:].strip().split(None,1)
                    #Each feature key could be used more than once,
                    #so store the entries as a list of strings.
                    if feature not in gf :
                        gf[feature] = [text]
                    else :
                        gf[feature].append(text)
                elif line[:5] == '#=GC ':
                    #Generic per-Column annotation, exactly 1 char per column
                    #Format: "#=GC <feature> <exactly 1 char per column>"
                    pass
                elif line[:5] == '#=GS ':
                    #Generic per-Sequence annotation, free text
                    #Format: "#=GS <seqname> <feature> <free text>"
                    id, feature, text = line[5:].strip().split(None,2)
                    #if id not in ids :
                    #    ids.append(id)
                    if id not in gs :
                        gs[id] = {}
                    if feature not in gs[id] :
                        gs[id][feature] = [text]
                    else :
                        gs[id][feature].append(text)
                elif line[:5] == "#=GR " :
                    #Generic per-Sequence AND per-Column markup
                    #Format: "#=GR <seqname> <feature> <exactly 1 char per column>"
                    id, feature, text = line[5:].strip().split(None,2)
                    #if id not in ids :
                    #    ids.append(id)
                    if id not in gr :
                        gr[id] = {}
                    if feature not in gr[id]:
                        gr[id][feature] = ""
                    gr[id][feature] += text.strip() # append to any previous entry
                    #TODO - Should we check the length matches the alignment length?
                    #       For iterlaced sequences the GR data can be split over
                    #       multiple lines
            #Next line...            


        assert len(seqs) <= len(ids)
        #assert len(gs)   <= len(ids)
        #assert len(gr)   <= len(ids)

        self.ids = ids
        self.sequences = seqs
        self.seq_annotation = gs
        self.seq_col_annotation = gr

        if ids and seqs :

            if self.records_per_alignment is not None \
            and self.records_per_alignment <> len(ids) :
                raise ValueError("Found %i records in this alignment, told to expect %i" \
                                 % (len(ids), self.records_per_alignment))

            alignment = Alignment(self.alphabet)

            #TODO - Introduce an annotated alignment class?
            #For now, store the annotation a new private property:
            alignment._annotations = gr

            alignment_length = len(seqs.values()[0])
            for id in ids :
                seq = seqs[id]
                if alignment_length <> len(seq) :
                    raise ValueError("Sequences have different lengths, or repeated identifier")
                name, start, end = self._identifier_split(id)
                alignment.add_sequence(id, seq, start=start, end=end)

                record = alignment.get_all_seqs()[-1]

                assert record.id == id or record.description == id
                
                record.id = id
                record.name = name
                record.description = id
                
                #will be overridden by _populate_meta_data if an explicit
                #accession is provided:
                record.annotations["accession"]=name

                self._populate_meta_data(id, record)
            return alignment
        else :
            return None


    def _identifier_split(self, identifier) :
        """Returns (name,start,end) string tuple from an identier."""
        if identifier.find("/")<>-1 :
            start_end = identifier.split("/",1)[1]
            if start_end.count("-")==1 :
                start, end = map(int, start_end.split("-"))
                name = identifier.split("/",1)[0]
                return (name, start, end)
        return (identifier, None, None)
    
    def _get_meta_data(self, identifier, meta_dict) :
        """Takes an itentifier and returns dict of all meta-data matching it.

        For example, given "Q9PN73_CAMJE/149-220" will return all matches to
        this or "Q9PN73_CAMJE" which the identifier without its /start-end
        suffix.

        In the example below, the suffix is required to match the AC, but must
        be removed to match the OS and OC meta-data.

        # STOCKHOLM 1.0
        #=GS Q9PN73_CAMJE/149-220  AC Q9PN73
        ...
        Q9PN73_CAMJE/149-220               NKA...
        ...
        #=GS Q9PN73_CAMJE OS Campylobacter jejuni
        #=GS Q9PN73_CAMJE OC Bacteria 

        This function will return an empty dictionary if no data is found."""
        name, start, end = self._identifier_split(identifier)
        if name==identifier :
            identifier_keys = [identifier]
        else :
            identifier_keys = [identifier, name]
        answer = {}
        for identifier_key in identifier_keys :
            try :
                for feature_key in meta_dict[identifier_key] :
                    answer[feature_key] = meta_dict[identifier_key][feature_key]
            except KeyError :
                pass
        return answer

    def _populate_meta_data(self, identifier, record) :
        """Adds meta-date to a SecRecord's annotations dictionary.

        This function applies the PFAM conventions."""

        seq_data = self._get_meta_data(identifier, self.seq_annotation)
        for feature in seq_data :
            #Note this dictionary contains lists!
            if feature=="AC" : #ACcession number
                assert len(seq_data[feature])==1
                record.annotations["accession"]=seq_data[feature][0]
            elif feature=="DE" : #DEscription
                record.description = "\n".join(seq_data[feature])
            elif feature=="DR" : #Database Reference
                #Should we try and parse the strings?
                record.dbxrefs = seq_data[feature]
            elif feature in self.pfam_gs_mapping :
                record.annotations[self.pfam_gs_mapping[feature]] = ", ".join(seq_data[feature])
            else :
                #Ignore it?
                record.annotations["GS:" + feature] = ", ".join(seq_data[feature])

        seq_col_data = self._get_meta_data(identifier, self.seq_col_annotation)
        for feature in seq_col_data :
            #Note this dictionary contains strings!
            if feature in self.pfam_gr_mapping :
                record.annotations[self.pfam_gr_mapping[feature]] = seq_col_data[feature]
            else :
                #Ignore it?
                record.annotations["GR:" + feature] = seq_col_data[feature]
    
if __name__ == "__main__" :
    print "Testing..."
    from cStringIO import StringIO

    # This example with its slightly odd (partial) annotation is from here:
    # http://www.cgb.ki.se/cgb/groups/sonnhammer/Stockholm.html
    # I don't know what the "GR_O31699/88-139_IN ..." line is meant to be.
    sth_example = \
"""# STOCKHOLM 1.0
#=GF ID CBS
#=GF AC PF00571
#=GF DE CBS domain
#=GF AU Bateman A
#=GF CC CBS domains are small intracellular modules mostly found  
#=GF CC in 2 or four copies within a protein. 
#=GF SQ 67
#=GS O31698/18-71 AC O31698
#=GS O83071/192-246 AC O83071
#=GS O83071/259-312 AC O83071
#=GS O31698/88-139 AC O31698
#=GS O31698/88-139 OS Bacillus subtilis
O83071/192-246          MTCRAQLIAVPRASSLAE..AIACAQKM....RVSRVPVYERS
#=GR O83071/192-246 SA  999887756453524252..55152525....36463774777
O83071/259-312          MQHVSAPVFVFECTRLAY..VQHKLRAH....SRAVAIVLDEY
#=GR O83071/259-312 SS  CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEE
O31698/18-71            MIEADKVAHVQVGNNLEH..ALLVLTKT....GYTAIPVLDPS
#=GR O31698/18-71 SS    CCCHHHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEHHH
O31698/88-139           EVMLTDIPRLHINDPIMK..GFGMVINN......GFVCVENDE
#=GR O31698/88-139 SS   CCCCCCCHHHHHHHHHHH..HEEEEEEE....EEEEEEEEEEH
#=GC SS_cons            CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEH
O31699/88-139           EVMLTDIPRLHINDPIMK..GFGMVINN......GFVCVENDE
#=GR O31699/88-139 AS   ________________*__________________________
#=GR_O31699/88-139_IN   ____________1______________2__________0____
//
"""

    # Interlaced example from BioPerl documentation.  Also note the blank line.
    # http://www.bioperl.org/wiki/Stockholm_multiple_alignment_format
    sth_example2 = \
"""# STOCKHOLM 1.0
#=GC SS_cons       .................<<<<<<<<...<<<<<<<........>>>>>>>..
AP001509.1         UUAAUCGAGCUCAACACUCUUCGUAUAUCCUC-UCAAUAUGG-GAUGAGGGU
#=GR AP001509.1 SS -----------------<<<<<<<<---..<<-<<-------->>->>..--
AE007476.1         AAAAUUGAAUAUCGUUUUACUUGUUUAU-GUCGUGAAU-UGG-CACGA-CGU
#=GR AE007476.1 SS -----------------<<<<<<<<-----<<.<<-------->>.>>----

#=GC SS_cons       ......<<<<<<<.......>>>>>>>..>>>>>>>>...............
AP001509.1         CUCUAC-AGGUA-CCGUAAA-UACCUAGCUACGAAAAGAAUGCAGUUAAUGU
#=GR AP001509.1 SS -------<<<<<--------->>>>>--->>>>>>>>---------------
AE007476.1         UUCUACAAGGUG-CCGG-AA-CACCUAACAAUAAGUAAGUCAGCAGUGAGAU
#=GR AE007476.1 SS ------.<<<<<--------->>>>>.-->>>>>>>>---------------
//"""

    print "--------"
    print "StockholmIterator(stockholm alignment file)"
    iterator = StockholmIterator(StringIO(sth_example))
    count=0
    for alignment in iterator :
        for record in alignment.get_all_seqs() :
            count=count+1
    assert count == 5
    
    #Check the last record...
    assert record.id == 'O31699/88-139'
    assert record.name == 'O31699'
    assert record.description == 'O31699/88-139'
    assert len(record.annotations)==4+1 #weight
    assert record.annotations["accession"]=='O31699'
    assert record.annotations["start"]==88
    assert record.annotations["end"]==139
    assert record.annotations["active_site"]=='________________*__________________________'

    iterator = StockholmIterator(StringIO(sth_example))
    count=0
    for alignment in iterator :
        for record in alignment.get_all_seqs() :
            count=count+1
            break
        break
    assert count==1
    #Check the last record...
    assert record.id == 'O83071/192-246'
    assert record.name == 'O83071'
    assert record.description == 'O83071/192-246'
    assert len(record.annotations)==4 + 1#weight
    assert record.annotations["accession"]=='O83071'
    assert record.annotations["start"]==192
    assert record.annotations["end"]==246
    assert record.annotations["surface_accessibility"]=="999887756453524252..55152525....36463774777"

    assert [[r.id for r in a.get_all_seqs()] \
            for a in StockholmIterator(StringIO(sth_example))] \
        == [['O83071/192-246', 'O83071/259-312', 'O31698/18-71', 'O31698/88-139', 'O31699/88-139']]

    assert [[r.name for r in a.get_all_seqs()] \
            for a in StockholmIterator(StringIO(sth_example))] \
        == [['O83071', 'O83071', 'O31698', 'O31698', 'O31699']]

    assert [[r.description for r in a.get_all_seqs()] \
            for a in StockholmIterator(StringIO(sth_example))] \
        == [['O83071/192-246', 'O83071/259-312', 'O31698/18-71', 'O31698/88-139', 'O31699/88-139']]


    print "--------"
    print "StockholmIterator(interlaced stockholm alignment file)"
    iterator = iter(StockholmIterator(StringIO(sth_example2)).next().get_all_seqs())
    record = iterator.next()
    assert record.id == "AP001509.1"
    assert len(record.seq) == 104
    assert "secondary_structure" in record.annotations
    assert len(record.annotations["secondary_structure"]) == 104
    record = iterator.next()
    assert record.id == "AE007476.1"
    assert len(record.seq) == 104
    assert "secondary_structure" in record.annotations
    assert len(record.annotations["secondary_structure"]) == 104
    try :
        record = iterator.next()
    except StopIteration :
        record = None
    assert record is None
    
    print "--------"
    print "StockholmIterator(concatenated stockholm alignment files)"
    assert 2 == len(list(StockholmIterator(StringIO(sth_example + sth_example2))))
    assert 2 == len(list(StockholmIterator(StringIO(sth_example + "\n" + sth_example2))))
    assert 2 == len(list(StockholmIterator(StringIO(sth_example + "\n\n" + sth_example2))))
    assert 2 == len(list(StockholmIterator(StringIO(sth_example2 + "\n\n" + sth_example))))
    assert 2 == len(list(StockholmIterator(StringIO(sth_example2 + "\n" + sth_example))))

    print "--------"
    print "Checking write/read"
    handle = StringIO()
    list1 = list(StockholmIterator(StringIO(sth_example2 + "\n" + sth_example)))
    StockholmWriter(handle).write_file(list1)
    handle.seek(0)
    list2 = list(StockholmIterator(handle))
    assert len(list1) == len(list2)
    for a1,a2 in zip(list1, list2) :
        assert len(a1.get_all_seqs()) == len(a2.get_all_seqs())
        for r1, r2 in zip(a1.get_all_seqs(), a2.get_all_seqs()) :
            assert r1.id == r2.id
            assert r1.seq.tostring() == r2.seq.tostring()
    print "Done"
