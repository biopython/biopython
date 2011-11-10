# Copyright 2008-2011 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Bio.AlignIO support for "fasta-m10" output from Bill Pearson's FASTA tools.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).

This module contains a parser for the pairwise alignments produced by Bill
Pearson's FASTA tools, for use from the Bio.AlignIO interface where it is
refered to as the "fasta-m10" file format (as we only support the machine
readable output format selected with the -m 10 command line option).

This module does NOT cover the generic "fasta" file format originally
developed as an input format to the FASTA tools.  The Bio.AlignIO and
Bio.SeqIO both use the Bio.SeqIO.FastaIO module to deal with these files,
which can also be used to store a multiple sequence alignments.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Interfaces import AlignmentIterator
from Bio.Alphabet import single_letter_alphabet, generic_dna, generic_protein
from Bio.Alphabet import Gapped


def _extract_alignment_region(alignment_seq_with_flanking, annotation):
    """Helper function for the main parsing code (PRIVATE).

    To get the actual pairwise alignment sequences, we must first
    translate the un-gapped sequence based coordinates into positions
    in the gapped sequence (which may have a flanking region shown
    using leading - characters).  To date, I have never seen any
    trailing flanking region shown in the m10 file, but the
    following code should also cope with that.

    Note that this code seems to work fine even when the "sq_offset"
    entries are prsent as a result of using the -X command line option.
    """
    align_stripped = alignment_seq_with_flanking.strip("-")
    display_start = int(annotation['al_display_start'])
    if int(annotation['al_start']) <= int(annotation['al_stop']):
        start = int(annotation['al_start']) \
              - display_start
        end   = int(annotation['al_stop']) \
              - display_start + 1
    else:
        #FASTA has flipped this sequence...
        start = display_start \
              - int(annotation['al_start'])
        end   = display_start \
              - int(annotation['al_stop']) + 1
    end += align_stripped.count("-")
    assert 0 <= start and start < end and end <= len(align_stripped), \
           "Problem with sequence start/stop,\n%s[%i:%i]\n%s" \
           % (alignment_seq_with_flanking, start, end, annotation)
    return align_stripped[start:end]

def FastaM10Iterator(handle, alphabet = single_letter_alphabet):
    """Alignment iterator for the FASTA tool's pairwise alignment output.

    This is for reading the pairwise alignments output by Bill Pearson's
    FASTA program when called with the -m 10 command line option for machine
    readable output.  For more details about the FASTA tools, see the website
    http://fasta.bioch.virginia.edu/ and the paper:

         W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448

    This class is intended to be used via the Bio.AlignIO.parse() function
    by specifying the format as "fasta-m10" as shown in the following code:

        from Bio import AlignIO
        handle = ...
        for a in AlignIO.parse(handle, "fasta-m10"):
            assert len(a) == 2, "Should be pairwise!"
            print "Alignment length %i" % a.get_alignment_length()
            for record in a:
                print record.seq, record.name, record.id

    Note that this is not a full blown parser for all the information
    in the FASTA output - for example, most of the header and all of the
    footer is ignored.  Also, the alignments are not batched according to
    the input queries.

    Also note that there can be up to about 30 letters of flanking region
    included in the raw FASTA output as contextual information.  This is NOT
    part of the alignment itself, and is not included in the resulting
    MultipleSeqAlignment objects returned.
    """
    if alphabet is None:
        alphabet = single_letter_alphabet
    
    state_PREAMBLE = -1
    state_NONE = 0
    state_QUERY_HEADER = 1
    state_ALIGN_HEADER = 2
    state_ALIGN_QUERY = 3
    state_ALIGN_MATCH = 4
    state_ALIGN_CONS = 5

    def build_hsp():
        if not query_tags and not match_tags:
            raise ValueError("No data for query %r, match %r" \
                             % (query_id, match_id))
        assert query_tags, query_tags
        assert match_tags, match_tags
        evalue = align_tags.get("fa_expect", None)
        q = "?" #Just for printing len(q) in debug below
        m = "?" #Just for printing len(m) in debug below
        tool = global_tags.get("tool", "").upper()
        try:
            q = _extract_alignment_region(query_seq, query_tags)
            if tool in ["TFASTX"] and len(match_seq) == len(q):
                m = match_seq
                #Quick hack until I can work out how -, * and / characters
                #and the apparent mix of aa and bp coordindates works.
            else:
                m = _extract_alignment_region(match_seq, match_tags)
            assert len(q) == len(m)
        except AssertionError, err:
            print "Darn... amino acids vs nucleotide coordinates?"
            print tool
            print query_seq
            print query_tags
            print q, len(q)
            print match_seq
            print match_tags
            print m, len(m)
            print handle.name
            raise err

        assert alphabet is not None
        alignment = MultipleSeqAlignment([], alphabet)

        #TODO - Introduce an annotated alignment class?
        #For now, store the annotation a new private property:
        alignment._annotations = {}
        
        #Want to record both the query header tags, and the alignment tags.
        for key, value in header_tags.iteritems():
            alignment._annotations[key] = value
        for key, value in align_tags.iteritems():
            alignment._annotations[key] = value
        
        #Query
        #=====
        record = SeqRecord(Seq(q, alphabet),
                           id = query_id,
                           name = "query",
                           description = query_descr,
                           annotations = {"original_length" : int(query_tags["sq_len"])})
        #TODO - handle start/end coordinates properly. Short term hack for now:
        record._al_start = int(query_tags["al_start"])
        record._al_stop = int(query_tags["al_stop"])
        alignment.append(record)

        #TODO - What if a specific alphabet has been requested?
        #TODO - Use an IUPAC alphabet?
        #TODO - Can FASTA output RNA?
        if alphabet == single_letter_alphabet and "sq_type" in query_tags:
            if query_tags["sq_type"] == "D":
                record.seq.alphabet = generic_dna
            elif query_tags["sq_type"] == "p":
                record.seq.alphabet = generic_protein
        if "-" in q:
            if not hasattr(record.seq.alphabet,"gap_char"):
                record.seq.alphabet = Gapped(record.seq.alphabet, "-")

        #Match
        #=====
        record = SeqRecord(Seq(m, alphabet),
                           id = match_id,
                           name = "match",
                           description = match_descr,
                           annotations = {"original_length" : int(match_tags["sq_len"])})
        #TODO - handle start/end coordinates properly. Short term hack for now:
        record._al_start = int(match_tags["al_start"])
        record._al_stop = int(match_tags["al_stop"])
        alignment.append(record)

        #This is still a very crude way of dealing with the alphabet:
        if alphabet == single_letter_alphabet and "sq_type" in match_tags:
            if match_tags["sq_type"] == "D":
                record.seq.alphabet = generic_dna
            elif match_tags["sq_type"] == "p":
                record.seq.alphabet = generic_protein
        if "-" in m:
            if not hasattr(record.seq.alphabet,"gap_char"):
                record.seq.alphabet = Gapped(record.seq.alphabet, "-")

        return alignment

    state = state_PREAMBLE
    query_id = None
    match_id = None
    query_descr = ""
    match_descr = ""
    global_tags = {}
    header_tags = {}
    align_tags = {}
    query_tags = {}
    match_tags = {}
    query_seq = ""
    match_seq = ""
    cons_seq = ""
    for line in handle:
        if ">>>" in line and not line.startswith(">>>"):
            if query_id and match_id:
                #This happens on old FASTA output which lacked an end of
                #query >>><<< marker line.
                yield build_hsp()
            state = state_NONE
            query_descr = line[line.find(">>>")+3:].strip()
            query_id = query_descr.split(None,1)[0]
            match_id = None
            header_tags = {}
            align_tags = {}
            query_tags = {}
            match_tags = {}
            query_seq = ""
            match_seq = ""
            cons_seq = ""
        elif line.startswith("!! No "):
            #e.g.
            #!! No library sequences with E() < 0.5
            #or on more recent versions,
            #No sequences with E() < 0.05
            assert state == state_NONE
            assert not header_tags
            assert not align_tags
            assert not match_tags
            assert not query_tags
            assert match_id is None
            assert not query_seq
            assert not match_seq
            assert not cons_seq
            query_id = None
        elif line.strip() in [">>><<<", ">>>///"]:
            #End of query, possible end of all queries
            if query_id and match_id:
                yield build_hsp()
            state = state_NONE
            query_id = None
            match_id = None
            header_tags = {}
            align_tags = {}
            query_tags = {}
            match_tags = {}
            query_seq = ""
            match_seq = ""
            cons_seq = ""
        elif line.startswith(">>>"):
            #Should be start of a match!
            assert query_id is not None
            assert line[3:].split(", ",1)[0] == query_id, line
            assert match_id is None
            assert not header_tags
            assert not align_tags
            assert not query_tags
            assert not match_tags
            assert not match_seq
            assert not query_seq
            assert not cons_seq
            state = state_QUERY_HEADER
        elif line.startswith(">>"):
            #Should now be at start of a match alignment!
            if query_id and match_id:
                yield build_hsp()
            align_tags = {}
            query_tags = {}
            match_tags = {}
            query_seq = ""
            match_seq = ""
            cons_seq = ""
            match_descr = line[2:].strip()
            match_id = match_descr.split(None,1)[0]
            state = state_ALIGN_HEADER
        elif line.startswith(">--"):
            #End of one HSP
            assert query_id and match_id, line
            yield build_hsp()
            #Clean up read for next HSP
            #but reuse header_tags
            align_tags = {}
            query_tags = {}
            match_tags = {}
            query_seq = ""
            match_seq = ""
            cons_seq = ""
            state = state_ALIGN_HEADER
        elif line.startswith(">"):
            if state == state_ALIGN_HEADER:
                #Should be start of query alignment seq...
                assert query_id is not None, line
                assert match_id is not None, line
                assert query_id.startswith(line[1:].split(None,1)[0]), line
                state = state_ALIGN_QUERY
            elif state == state_ALIGN_QUERY:
                #Should be start of match alignment seq
                assert query_id is not None, line
                assert match_id is not None, line
                assert match_id.startswith(line[1:].split(None,1)[0]), line
                state = state_ALIGN_MATCH
            elif state == state_NONE:
                #Can get > as the last line of a histogram
                pass
            else:
                assert False, "state %i got %r" % (state, line)
        elif line.startswith("; al_cons"):
            assert state == state_ALIGN_MATCH, line
            state = state_ALIGN_CONS
            #Next line(s) should be consensus seq...
        elif line.startswith("; "):
            if ": " in line:
                key, value = [s.strip() for s in line[2:].split(": ",1)]
            else:
                import warnings
                #Seen in lalign36, specifically version 36.3.4 Apr, 2011
                #Fixed in version 36.3.5b Oct, 2011(preload8)
                warnings.warn("Missing colon in line: %r" % line)
                try:
                    key, value = [s.strip() for s in line[2:].split(" ",1)]
                except ValueError:
                    raise ValueError("Bad line: %r" % line)
            if state == state_QUERY_HEADER:
                header_tags[key] = value
            elif state == state_ALIGN_HEADER:
                align_tags[key] = value
            elif state == state_ALIGN_QUERY:
                query_tags[key] = value
            elif state == state_ALIGN_MATCH:
                match_tags[key] = value
            else:
                assert False, "Unexpected state %r, %r" % (state, line)
        elif state == state_ALIGN_QUERY:
            query_seq += line.strip()
        elif state == state_ALIGN_MATCH:
            match_seq += line.strip()
        elif state == state_ALIGN_CONS:
            cons_seq += line.strip("\n")
        elif state == state_PREAMBLE:
            if line.startswith("#"):
                global_tags["command"] = line[1:].strip()
            elif line.startswith(" version "):
                global_tags["version"] = line[9:].strip()
            elif " compares a " in line:
                global_tags["tool"] = line[:line.find(" compares a ")].strip()
            elif " searches a " in line:
                global_tags["tool"] = line[:line.find(" searches a ")].strip()
        else:
            pass


if __name__ == "__main__":
    print "Running a quick self-test"

    #http://emboss.sourceforge.net/docs/themes/alnformats/align.simple
    simple_example = \
"""# /opt/fasta/fasta34 -Q -H -E 1 -m 10 NC_002127.faa NC_009649.faa
FASTA searches a protein or DNA sequence data bank
 version 34.26 January 12, 2007
Please cite:
 W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448

Query library NC_002127.faa vs NC_009649.faa library
searching NC_009649.faa library

  1>>>gi|10955263|ref|NP_052604.1| plasmid mobilization [Escherichia coli O157:H7 s 107 aa - 107 aa
 vs  NC_009649.faa library

  45119 residues in   180 sequences
  Expectation_n fit: rho(ln(x))= 6.9146+/-0.0249; mu= -5.7948+/- 1.273
 mean_var=53.6859+/-13.609, 0's: 0 Z-trim: 1  B-trim: 9 in 1/25
 Lambda= 0.175043

FASTA (3.5 Sept 2006) function [optimized, BL50 matrix (15:-5)] ktup: 2
 join: 36, opt: 24, open/ext: -10/-2, width:  16
 Scan time:  0.000
The best scores are:                                      opt bits E(180)
gi|152973457|ref|YP_001338508.1| ATPase with chape ( 931)   71 24.9    0.58
gi|152973588|ref|YP_001338639.1| F pilus assembly  ( 459)   63 23.1    0.99

>>>gi|10955263|ref|NP_052604.1|, 107 aa vs NC_009649.faa library
; pg_name: /opt/fasta/fasta34
; pg_ver: 34.26
; pg_argv: /opt/fasta/fasta34 -Q -H -E 1 -m 10 NC_002127.faa NC_009649.faa
; pg_name: FASTA
; pg_ver: 3.5 Sept 2006
; pg_matrix: BL50 (15:-5)
; pg_open-ext: -10 -2
; pg_ktup: 2
; pg_optcut: 24
; pg_cgap: 36
; mp_extrap: 60000 180
; mp_stats:  Expectation_n fit: rho(ln(x))= 6.9146+/-0.0249; mu= -5.7948+/- 1.273  mean_var=53.6859+/-13.609, 0's: 0 Z-trim: 1  B-trim: 9 in 1/25  Lambda= 0.175043
; mp_KS: -0.0000 (N=0) at 8159228
>>gi|152973457|ref|YP_001338508.1| ATPase with chaperone activity, ATP-binding subunit [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]
; fa_frame: f
; fa_initn:  65
; fa_init1:  43
; fa_opt:  71
; fa_z-score: 90.3
; fa_bits: 24.9
; fa_expect:   0.58
; sw_score: 71
; sw_ident: 0.250
; sw_sim: 0.574
; sw_overlap: 108
>gi|10955263| ..
; sq_len: 107
; sq_offset: 1
; sq_type: p
; al_start: 5
; al_stop: 103
; al_display_start: 1
--------------------------MTKRSGSNT-RRRAISRPVRLTAE
ED---QEIRKRAAECGKTVSGFLRAAALGKKVNSLTDDRVLKEVM-----
RLGALQKKLFIDGKRVGDREYAEVLIAITEYHRALLSRLMAD
>gi|152973457|ref|YP_001338508.1| ..
; sq_len: 931
; sq_type: p
; al_start: 96
; al_stop: 195
; al_display_start: 66
SDFFRIGDDATPVAADTDDVVDASFGEPAAAGSGAPRRRGSGLASRISEQ
SEALLQEAAKHAAEFGRS------EVDTEHLLLALADSDVVKTILGQFKI
KVDDLKRQIESEAKR-GDKPF-EGEIGVSPRVKDALSRAFVASNELGHSY
VGPEHFLIGLAEEGEGLAANLLRRYGLTPQ
>>gi|152973588|ref|YP_001338639.1| F pilus assembly protein [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]
; fa_frame: f
; fa_initn:  33
; fa_init1:  33
; fa_opt:  63
; fa_z-score: 86.1
; fa_bits: 23.1
; fa_expect:   0.99
; sw_score: 63
; sw_ident: 0.266
; sw_sim: 0.656
; sw_overlap: 64
>gi|10955263| ..
; sq_len: 107
; sq_offset: 1
; sq_type: p
; al_start: 32
; al_stop: 94
; al_display_start: 2
TKRSGSNTRRRAISRPVRLTAEEDQEIRKRAAECGKTVSGFLRAAALGKK
VNSLTDDRVLKEV-MRLGALQKKLFIDGKRVGDREYAEVLIAITEYHRAL
LSRLMAD
>gi|152973588|ref|YP_001338639.1| ..
; sq_len: 459
; sq_type: p
; al_start: 191
; al_stop: 248
; al_display_start: 161
VGGLFPRTQVAQQKVCQDIAGESNIFSDWAASRQGCTVGG--KMDSVQDK
ASDKDKERVMKNINIMWNALSKNRLFDG----NKELKEFIMTLTGTLIFG
ENSEITPLPARTTDQDLIRAMMEGGTAKIYHCNDSDKCLKVVADATVTIT
SNKALKSQISALLSSIQNKAVADEKLTDQE
  2>>>gi|10955264|ref|NP_052605.1| hypothetical protein pOSAK1_02 [Escherichia coli O157:H7 s 126 aa - 126 aa
 vs  NC_009649.faa library

  45119 residues in   180 sequences
  Expectation_n fit: rho(ln(x))= 7.1374+/-0.0246; mu= -7.6540+/- 1.313
 mean_var=51.1189+/-13.171, 0's: 0 Z-trim: 1  B-trim: 8 in 1/25
 Lambda= 0.179384

FASTA (3.5 Sept 2006) function [optimized, BL50 matrix (15:-5)] ktup: 2
 join: 36, opt: 24, open/ext: -10/-2, width:  16
 Scan time:  0.000
The best scores are:                                      opt bits E(180)
gi|152973462|ref|YP_001338513.1| hypothetical prot ( 101)   58 22.9    0.29

>>>gi|10955264|ref|NP_052605.1|, 126 aa vs NC_009649.faa library
; pg_name: /opt/fasta/fasta34
; pg_ver: 34.26
; pg_argv: /opt/fasta/fasta34 -Q -H -E 1 -m 10 NC_002127.faa NC_009649.faa
; pg_name: FASTA
; pg_ver: 3.5 Sept 2006
; pg_matrix: BL50 (15:-5)
; pg_open-ext: -10 -2
; pg_ktup: 2
; pg_optcut: 24
; pg_cgap: 36
; mp_extrap: 60000 180
; mp_stats:  Expectation_n fit: rho(ln(x))= 7.1374+/-0.0246; mu= -7.6540+/- 1.313  mean_var=51.1189+/-13.171, 0's: 0 Z-trim: 1  B-trim: 8 in 1/25  Lambda= 0.179384
; mp_KS: -0.0000 (N=0) at 8159228
>>gi|152973462|ref|YP_001338513.1| hypothetical protein KPN_pKPN3p05904 [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]
; fa_frame: f
; fa_initn:  50
; fa_init1:  50
; fa_opt:  58
; fa_z-score: 95.8
; fa_bits: 22.9
; fa_expect:   0.29
; sw_score: 58
; sw_ident: 0.289
; sw_sim: 0.632
; sw_overlap: 38
>gi|10955264| ..
; sq_len: 126
; sq_offset: 1
; sq_type: p
; al_start: 1
; al_stop: 38
; al_display_start: 1
------------------------------MKKDKKYQIEAIKNKDKTLF
IVYATDIYSPSEFFSKIESDLKKKKSKGDVFFDLIIPNGGKKDRYVYTSF
NGEKFSSYTLNKVTKTDEYN
>gi|152973462|ref|YP_001338513.1| ..
; sq_len: 101
; sq_type: p
; al_start: 44
; al_stop: 81
; al_display_start: 14
DALLGEIQRLRKQVHQLQLERDILTKANELIKKDLGVSFLKLKNREKTLI
VDALKKKYPVAELLSVLQLARSCYFYQNVCTISMRKYA
  3>>>gi|10955265|ref|NP_052606.1| hypothetical protein pOSAK1_03 [Escherichia coli O157:H7 s 346 aa - 346 aa
 vs  NC_009649.faa library

  45119 residues in   180 sequences
  Expectation_n fit: rho(ln(x))= 6.0276+/-0.0276; mu= 3.0670+/- 1.461
 mean_var=37.1634+/- 8.980, 0's: 0 Z-trim: 1  B-trim: 14 in 1/25
 Lambda= 0.210386

FASTA (3.5 Sept 2006) function [optimized, BL50 matrix (15:-5)] ktup: 2
 join: 37, opt: 25, open/ext: -10/-2, width:  16
 Scan time:  0.020
The best scores are:                                      opt bits E(180)
gi|152973545|ref|YP_001338596.1| putative plasmid  ( 242)   70 27.5   0.082

>>>gi|10955265|ref|NP_052606.1|, 346 aa vs NC_009649.faa library
; pg_name: /opt/fasta/fasta34
; pg_ver: 34.26
; pg_argv: /opt/fasta/fasta34 -Q -H -E 1 -m 10 NC_002127.faa NC_009649.faa
; pg_name: FASTA
; pg_ver: 3.5 Sept 2006
; pg_matrix: BL50 (15:-5)
; pg_open-ext: -10 -2
; pg_ktup: 2
; pg_optcut: 25
; pg_cgap: 37
; mp_extrap: 60000 180
; mp_stats:  Expectation_n fit: rho(ln(x))= 6.0276+/-0.0276; mu= 3.0670+/- 1.461  mean_var=37.1634+/- 8.980, 0's: 0 Z-trim: 1  B-trim: 14 in 1/25  Lambda= 0.210386
; mp_KS: -0.0000 (N=0) at 8159228
>>gi|152973545|ref|YP_001338596.1| putative plasmid SOS inhibition protein A [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]
; fa_frame: f
; fa_initn:  52
; fa_init1:  52
; fa_opt:  70
; fa_z-score: 105.5
; fa_bits: 27.5
; fa_expect:  0.082
; sw_score: 70
; sw_ident: 0.279
; sw_sim: 0.651
; sw_overlap: 43
>gi|10955265| ..
; sq_len: 346
; sq_offset: 1
; sq_type: p
; al_start: 197
; al_stop: 238
; al_display_start: 167
DFMCSILNMKEIVEQKNKEFNVDIKKETIESELHSKLPKSIDKIHEDIKK
QLSC-SLIMKKIDVEMEDYSTYCFSALRAIEGFIYQILNDVCNPSSSKNL
GEYFTENKPKYIIREIHQET
>gi|152973545|ref|YP_001338596.1| ..
; sq_len: 242
; sq_type: p
; al_start: 52
; al_stop: 94
; al_display_start: 22
IMTVEEARQRGARLPSMPHVRTFLRLLTGCSRINSDVARRIPGIHRDPKD
RLSSLKQVEEALDMLISSHGEYCPLPLTMDVQAENFPEVLHTRTVRRLKR
QDFAFTRKMRREARQVEQSW
>>><<<


579 residues in 3 query   sequences
45119 residues in 180 library sequences
 Scomplib [34.26]
 start: Tue May 20 16:38:45 2008 done: Tue May 20 16:38:45 2008
 Total Scan time:  0.020 Total Display time:  0.010

Function used was FASTA [version 34.26 January 12, 2007]

"""                 


    from StringIO import StringIO

    alignments = list(FastaM10Iterator(StringIO(simple_example)))
    assert len(alignments) == 4, len(alignments)
    assert len(alignments[0]) == 2
    for a in alignments:
        print "Alignment %i sequences of length %i" \
              % (len(a), a.get_alignment_length())
        for r in a:
            print "%s %s %i" % (r.seq, r.id, r.annotations["original_length"])
        #print a.annotations
    print "Done"

    import os
    path = "../../Tests/Fasta/"
    files = [f for f in os.listdir(path) if os.path.splitext(f)[-1] == ".m10"]
    files.sort()
    for filename in files:
        if os.path.splitext(filename)[-1] == ".m10":
            print
            print filename
            print "="*len(filename)
            for i,a in enumerate(FastaM10Iterator(open(os.path.join(path,filename)))):
                print "#%i, %s" % (i+1,a)
                for r in a:
                    if "-" in r.seq:
                        assert r.seq.alphabet.gap_char == "-"
                    else:
                        assert not hasattr(r.seq.alphabet, "gap_char")
