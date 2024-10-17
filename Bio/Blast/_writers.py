import codecs
import html
from abc import ABC
from abc import abstractmethod

from Bio.Seq import UndefinedSequenceError


def _html_entity_replace(error):
    start = error.start
    end = error.end
    character = error.object[start:end]
    codepoint = ord(character)
    entity = f"&amp;{html.entities.codepoint2name[codepoint]};"
    return entity, end


codecs.register_error("htmlentityreplace", _html_entity_replace)


class BaseXMLWriter(ABC):
    """XML Writer, abstract base class."""

    def __init__(self, stream):
        """Init the stuff."""
        self.stream = stream

    def write(self, records):
        """Write the records."""
        stream = self.stream
        program = records.program
        self._program = program
        self._write_xml_declaration()
        self._write_definition()
        self._start_blastoutput()
        self._write_program(program.encode())
        self._write_version(records.version.encode())
        reference = html.escape(records.reference).encode("ASCII", "htmlentityreplace")
        self._write_reference(reference)
        self._write_db(records.db.encode())
        try:
            query = records.query
        except AttributeError:  # XML2
            pass
        else:
            self._write_query_id(query.id.encode())
            self._write_query_def(query.description.encode())
            self._write_query_len(len(query))
            try:
                query_seq = bytes(query.seq)
            except UndefinedSequenceError:
                pass
            else:
                self._write_query_seq(query_seq)
        self._write_params(records.param)
        count = self._write_records(records)
        try:
            mbstat = records.mbstat
        except AttributeError:
            pass
        else:
            self._start_mbstat()
            self._write_statistics(mbstat)
            self._end_mbstat()
        self._end_blastoutput()
        return count

    def _write_xml_declaration(self):
        self.stream.write(b'<?xml version="1.0"?>\n')

    def _write_params(self, param):
        self._start_param()
        self.stream.write(b"    <Parameters>\n")
        value = param.get("matrix")
        if value is not None:
            self._write_parameters_matrix(value.encode())
        self._write_parameters_expect(param["expect"])
        value = param.get("include")
        if value is not None:
            self._write_parameters_include(value)
        value = param.get("sc-match")
        if value is not None:
            self._write_parameters_sc_match(value)
        value = param.get("sc-mismatch")
        if value is not None:
            self._write_parameters_sc_mismatch(value)
        value = param.get("gap-open")
        if value is not None:
            self._write_parameters_gap_open(value)
        value = param.get("gap-extend")
        if value is not None:
            self._write_parameters_gap_extend(value)
        value = param.get("filter")
        if value is not None:
            self._write_parameters_filter(value.encode())
        value = param.get("pattern")
        if value is not None:
            self._write_parameters_pattern(value.encode())
        value = param.get("entrez-query")
        if value is not None:
            self._write_parameters_entrez_query(value.encode())
        value = param.get("cbs")
        if value is not None:
            self._write_parameters_cbs(value)
        value = param.get("query-gencode")
        if value is not None:
            self._write_parameters_query_gencode(value)
        value = param.get("db-gencode")
        if value is not None:
            self._write_parameters_db_gencode(value)
        value = param.get("bl2seq-mode")
        if value is not None:
            self._write_parameters_bl2seq_mode(value.encode())
        self.stream.write(b"    </Parameters>\n")
        self._end_param()

    def _write_records(self, records):
        count = 0
        self._start_iterations()
        for record in records:
            self._write_record(record)
            count += 1
        self._end_iterations()
        return count

    def _write_record(self, record):
        stream = self.stream
        self._start_iteration()
        try:
            num = record.num
        except AttributeError:  # XML2
            pass
        else:
            self._write_iteration_num(num)
        query = record.query
        if query is None:
            query_length = None
        else:
            query_length = len(query.seq)
            self._write_iteration_query_id(query.id.encode())
            self._write_iteration_query_def(query.description.encode())
            self._write_iteration_query_len(query_length)
            for feature in query.features:
                if feature.type == "masking":
                    self._write_iteration_query_masking(feature.location)
        self._start_iteration_hits()
        for hit in record:
            self._write_hit(hit, query_length)
        self._end_iteration_hits()
        try:
            stat = record.stat
        except AttributeError:
            pass
        else:
            self._start_iteration_stat()
            self._write_statistics(stat)
            self._end_iteration_stat()
        self._end_iteration()

    def _write_hit(self, hit, query_length):
        stream = self.stream
        program = self._program
        target = hit.target
        target_length = len(target.seq)
        self._start_iteration_hit()
        self._write_hit_num(hit.num)
        try:
            targets = hit.targets
        except AttributeError:  # XML
            self._write_hit_id(target.id.encode())
            self._write_hit_def(target.description.encode())
            self._write_hit_accession(target.name.encode())
        else:  # XML2
            self._start_hit_targets()
            for target in targets:
                self._start_hitdescr()
                self._write_hit_id(target.id.encode())
                self._write_hit_accession(target.name.encode())
                self._write_hit_def(target.description.encode())
                taxid = target.annotations.get("taxid")
                if taxid is not None:
                    self._write_hit_taxid(taxid)
                sciname = target.annotations.get("sciname")
                if sciname is not None:
                    self._write_hit_sciname(sciname.encode())
                self._end_hitdescr()
            self._end_hit_targets()
        self._write_hit_len(target_length)
        self._start_hit_hsps()
        for hsp in hit:
            self._write_hsp(hsp, query_length, target_length)
        self._end_hit_hsps()
        self._end_iteration_hit()

    def _write_hsp(self, hsp, query_length, target_length):
        stream = self.stream
        program = self._program
        query = hsp.query
        target = hsp.target
        coordinates = hsp.coordinates
        hit_from, query_from = coordinates[:, 0]
        hit_to, query_to = coordinates[:, -1]
        if program in ("blastn", "megablast"):
            if hit_from <= hit_to:
                hit_frame = 1
                hit_from += 1
            else:
                hit_frame = -1
                hit_to += 1
        elif program in ("blastp", "blastx", "rpsblast"):
            hit_from += 1
            hit_frame = 0
        elif program in ("tblastn", "tblastx"):
            feature = target.features[0]
            coded_by = feature.qualifiers["coded_by"]
            if coded_by.startswith("complement("):
                assert coded_by.endswith(")")
                coded_by = coded_by[11:-1]
                strand = -1
            else:
                strand = +1
            hit_id, hit_from_to = coded_by.split(":")
            hit_from, hit_to = hit_from_to.split("..")
            hit_from = int(hit_from)
            hit_to = int(hit_to)
            hit_start = hit_from - 1
            hit_end = hit_to
            if strand == +1:
                hit_frame = hit_start % 3 + 1
            else:
                hit_frame = (hit_end - target_length) % -3 - 1
        if program in ("blastn", "megablast"):
            if query_from <= query_to:
                query_from += 1
                query_frame = 1
            else:
                query_to += 1
                query_frame = -1
        elif program in ("blastp", "tblastn", "rpsblast"):
            query_from += 1
            query_frame = 0
        elif program in ("blastx", "tblastx"):
            feature = query.features[0]
            coded_by = feature.qualifiers["coded_by"]
            if coded_by.startswith("complement("):
                assert coded_by.endswith(")")
                coded_by = coded_by[11:-1]
                strand = -1
            else:
                strand = +1
            query_id, query_from_to = coded_by.split(":")
            query_from, query_to = query_from_to.split("..")
            query_from = int(query_from)
            query_to = int(query_to)
            query_start = query_from - 1
            query_end = query_to
            if strand == +1:
                query_frame = query_start % 3 + 1
            else:
                query_frame = (query_end - query_length) % -3 - 1
        hseq = hsp[0]
        qseq = hsp[1]
        align_len = len(hseq)
        annotations = hsp.annotations
        bit_score = annotations["bit score"]
        evalue = annotations["evalue"]
        identity = annotations["identity"]
        positive = annotations.get("positive")
        gaps = annotations.get("gaps")
        midline = annotations["midline"]
        self._start_hsp()
        self._write_hsp_num(hsp.num)
        self._write_hsp_bit_score(str(bit_score).encode())
        self._write_hsp_score(hsp.score)
        self._write_hsp_evalue(str(evalue).encode())
        self._write_hsp_query_from(query_from)
        self._write_hsp_query_to(query_to)
        self._write_hsp_hit_from(hit_from)
        self._write_hsp_hit_to(hit_to)
        self._write_hsp_query_frame(query_frame)
        self._write_hsp_hit_frame(hit_frame)
        self._write_hsp_identity(identity)
        if positive is not None:
            self._write_hsp_positive(positive)
        if gaps is not None:
            self._write_hsp_gaps(gaps)
        self._write_hsp_align_len(align_len)
        self._write_hsp_qseq(qseq.encode())
        self._write_hsp_hseq(hseq.encode())
        self._write_hsp_midline(midline.encode())
        self._end_hsp()

    def _write_statistics(self, stat):
        self.stream.write(b"    <Statistics>\n")
        self._write_statistics_db_num(stat["db-num"])
        self._write_statistics_db_len(stat["db-len"])
        self._write_statistics_hsp_len(stat["hsp-len"])
        self._write_statistics_eff_space(str(stat["eff-space"]).encode())
        self._write_statistics_kappa(stat["kappa"])
        self._write_statistics_lambda(stat["lambda"])
        self._write_statistics_entropy(stat["entropy"])
        self.stream.write(b"    </Statistics>\n")

    def _write_definition(self):
        return

    @abstractmethod
    def _start_blastoutput(self):
        return

    @abstractmethod
    def _end_blastoutput(self):
        return

    @abstractmethod
    def _write_program(self, program):
        return

    @abstractmethod
    def _write_version(self, version):
        return

    @abstractmethod
    def _write_reference(self, reference):
        return

    @abstractmethod
    def _write_db(self, db):
        return

    @abstractmethod
    def _start_param(self):
        return

    @abstractmethod
    def _end_param(self):
        return

    @abstractmethod
    def _start_iterations(self):
        return

    @abstractmethod
    def _end_iterations(self):
        return

    @abstractmethod
    def _start_iteration(self):
        return

    @abstractmethod
    def _end_iteration(self):
        return

    @abstractmethod
    def _write_iteration_query_id(self, query_id):
        return

    @abstractmethod
    def _write_iteration_query_def(self, query_def):
        return

    @abstractmethod
    def _write_iteration_query_len(self, query_len):
        return

    def _write_iteration_query_masking(self, location):
        return

    @abstractmethod
    def _start_iteration_hits(self):
        return

    @abstractmethod
    def _end_iteration_hits(self):
        return

    @abstractmethod
    def _start_iteration_stat(self):
        return

    @abstractmethod
    def _end_iteration_stat(self):
        return

    @abstractmethod
    def _write_parameters_matrix(self, value):
        return

    @abstractmethod
    def _write_parameters_expect(self, value):
        return

    @abstractmethod
    def _write_parameters_include(self, value):
        return

    @abstractmethod
    def _write_parameters_sc_match(self, value):
        return

    @abstractmethod
    def _write_parameters_sc_mismatch(self, value):
        return

    @abstractmethod
    def _write_parameters_gap_open(self, value):
        return

    @abstractmethod
    def _write_parameters_gap_extend(self, value):
        return

    @abstractmethod
    def _write_parameters_filter(self, value):
        return

    @abstractmethod
    def _write_parameters_pattern(self, value):
        return

    @abstractmethod
    def _write_parameters_entrez_query(self, value):
        return

    def _write_parameters_cbs(self, value):
        return

    def _write_parameters_query_gencode(self, value):
        return

    def _write_parameters_db_gencode(self, value):
        return

    def _write_parameters_bl2seq_mode(self, value):
        return

    @abstractmethod
    def _write_statistics_db_num(self, db_num):
        return

    @abstractmethod
    def _write_statistics_db_len(self, db_len):
        return

    @abstractmethod
    def _write_statistics_hsp_len(self, hsp_len):
        return

    @abstractmethod
    def _write_statistics_eff_space(self, value):
        return

    @abstractmethod
    def _write_statistics_kappa(self, value):
        return

    @abstractmethod
    def _write_statistics_lambda(self, value):
        return

    @abstractmethod
    def _write_statistics_entropy(self, value):
        return

    @abstractmethod
    def _start_iteration_hit(self):
        return

    @abstractmethod
    def _end_iteration_hit(self):
        return

    @abstractmethod
    def _write_hit_num(self, num):
        return

    @abstractmethod
    def _write_hit_id(self, hit_id):
        return

    @abstractmethod
    def _write_hit_def(self, hit_def):
        return

    @abstractmethod
    def _write_hit_accession(self, hit_accession):
        return

    def _write_hit_taxid(self, taxid):
        return

    def _write_hit_sciname(self, sciname):
        return

    @abstractmethod
    def _write_hit_len(self, hit_length):
        return

    @abstractmethod
    def _start_hit_hsps(self):
        return

    @abstractmethod
    def _end_hit_hsps(self):
        return

    @abstractmethod
    def _start_hsp(self):
        return

    @abstractmethod
    def _end_hsp(self):
        return

    @abstractmethod
    def _write_hsp_num(self, num):
        return

    @abstractmethod
    def _write_hsp_bit_score(self, bit_score):
        return

    @abstractmethod
    def _write_hsp_score(self, score):
        return

    @abstractmethod
    def _write_hsp_evalue(self, evalue):
        return

    @abstractmethod
    def _write_hsp_query_from(self, query_from):
        return

    @abstractmethod
    def _write_hsp_query_to(self, query_to):
        return

    @abstractmethod
    def _write_hsp_hit_from(self, hit_from):
        return

    @abstractmethod
    def _write_hsp_hit_to(self, hit_to):
        return

    @abstractmethod
    def _write_hsp_pattern_from(self, hit_from):
        return

    @abstractmethod
    def _write_hsp_pattern_to(self, hit_to):
        return

    @abstractmethod
    def _write_hsp_query_frame(self, query_frame):
        return

    @abstractmethod
    def _write_hsp_hit_frame(self, hit_frame):
        return

    @abstractmethod
    def _write_hsp_identity(self, identity):
        return

    @abstractmethod
    def _write_hsp_positive(self, positive):
        return

    @abstractmethod
    def _write_hsp_gaps(self, gaps):
        return

    @abstractmethod
    def _write_hsp_align_len(self, align_len):
        return

    @abstractmethod
    def _write_hsp_density(self, density):
        return

    @abstractmethod
    def _write_hsp_qseq(self, qseq):
        return

    @abstractmethod
    def _write_hsp_hseq(self, hseq):
        return

    @abstractmethod
    def _write_hsp_midline(self, midline):
        return


class XMLWriter(BaseXMLWriter):
    """XML Writer."""

    def _write_definition(self):
        self.stream.write(
            b"""\
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
"""
        )

    def _start_blastoutput(self):
        self.stream.write(b"<BlastOutput>\n")

    def _end_blastoutput(self):
        self.stream.write(b"</BlastOutput>\n")

    def _write_program(self, program):
        self.stream.write(b"<BlastOutput_program>%b</BlastOutput_program>\n" % program)

    def _write_version(self, version):
        self.stream.write(b"<BlastOutput_version>%b</BlastOutput_version>\n" % version)

    def _write_reference(self, reference):
        self.stream.write(
            b"<BlastOutput_reference>%b</BlastOutput_reference>\n" % reference
        )

    def _write_db(self, db):
        self.stream.write(b"<BlastOutput_db>%b</BlastOutput_db>\n" % db)

    def _write_query_id(self, query_id):
        self.stream.write(
            b"<BlastOutput_query-ID>%b</BlastOutput_query-ID>\n" % query_id
        )

    def _write_query_def(self, query_def):
        self.stream.write(
            b"<BlastOutput_query-def>%b</BlastOutput_query-def>\n" % query_def
        )

    def _write_query_len(self, query_len):
        self.stream.write(
            b"<BlastOutput_query-len>%d</BlastOutput_query-len>\n" % query_len
        )

    def _write_query_seq(self, query_seq):
        self.stream.write(
            b"<BlastOutput_query-seq>%b</BlastOutput_query-seq>\n" % query_seq
        )

    def _start_param(self):
        self.stream.write(b"  <BlastOutput_param>\n")

    def _end_param(self):
        self.stream.write(b"  </BlastOutput_param>\n")

    def _start_iterations(self):
        self.stream.write(b"<BlastOutput_iterations>\n")

    def _end_iterations(self):
        self.stream.write(b"</BlastOutput_iterations>\n")

    def _start_mbstat(self):
        self.stream.write(b"  <BlastOutput_mbstat>\n")

    def _end_mbstat(self):
        self.stream.write(b"  </BlastOutput_mbstat>\n")

    def _start_iteration(self):
        self.stream.write(b"<Iteration>\n")

    def _end_iteration(self):
        self.stream.write(b"</Iteration>\n")

    def _write_iteration_num(self, num):
        self.stream.write(b"  <Iteration_iter-num>%d</Iteration_iter-num>\n" % num)

    def _write_iteration_query_id(self, query_id):
        self.stream.write(b"<Iteration_query-ID>%b</Iteration_query-ID>\n" % query_id)

    def _write_iteration_query_def(self, query_def):
        self.stream.write(
            b"<Iteration_query-def>%b</Iteration_query-def>\n" % query_def
        )

    def _write_iteration_query_len(self, query_len):
        self.stream.write(
            b"<Iteration_query-len>%d</Iteration_query-len>\n" % query_len
        )

    def _start_iteration_hits(self):
        self.stream.write(b"<Iteration_hits>\n")

    def _end_iteration_hits(self):
        self.stream.write(b"</Iteration_hits>\n")

    def _start_iteration_stat(self):
        self.stream.write(b"  <Iteration_stat>\n")

    def _end_iteration_stat(self):
        self.stream.write(b"  </Iteration_stat>\n")

    def _write_parameters_matrix(self, value):
        self.stream.write(b"      <Parameters_matrix>%b</Parameters_matrix>\n" % value)

    def _write_parameters_expect(self, value):
        self.stream.write(b"      <Parameters_expect>%g</Parameters_expect>\n" % value)

    def _write_parameters_include(self, value):
        self.stream.write(
            b"      <Parameters_include>%g</Parameters_include>\n" % value
        )

    def _write_parameters_sc_match(self, value):
        self.stream.write(
            b"       <Parameters_sc-match>%d</Parameters_sc-match>\n" % value
        )

    def _write_parameters_sc_mismatch(self, value):
        self.stream.write(
            b"       <Parameters_sc-mismatch>%d</Parameters_sc-mismatch>\n" % value
        )

    def _write_parameters_gap_open(self, value):
        self.stream.write(
            b"       <Parameters_gap-open>%d</Parameters_gap-open>\n" % value
        )

    def _write_parameters_gap_extend(self, value):
        self.stream.write(
            b"       <Parameters_gap-extend>%d</Parameters_gap-extend>\n" % value
        )

    def _write_parameters_filter(self, value):
        self.stream.write(b"       <Parameters_filter>%b</Parameters_filter>\n" % value)

    def _write_parameters_pattern(self, value):
        self.stream.write(
            b"       <Parameters_pattern>%b</Parameters_pattern>\n" % value
        )

    def _write_parameters_entrez_query(self, value):
        self.stream.write(
            b"       <Parameters_entrez-query>%b</Parameters_entrez-query>\n" % value
        )

    def _write_statistics_db_num(self, db_num):
        self.stream.write(b"      <Statistics_db-num>%d</Statistics_db-num>\n" % db_num)

    def _write_statistics_db_len(self, db_len):
        self.stream.write(b"      <Statistics_db-len>%d</Statistics_db-len>\n" % db_len)

    def _write_statistics_hsp_len(self, hsp_len):
        self.stream.write(
            b"      <Statistics_hsp-len>%d</Statistics_hsp-len>\n" % hsp_len
        )

    def _write_statistics_eff_space(self, value):
        self.stream.write(
            b"      <Statistics_eff-space>%s</Statistics_eff-space>\n" % value
        )

    def _write_statistics_kappa(self, value):
        self.stream.write(b"      <Statistics_kappa>%r</Statistics_kappa>\n" % value)

    def _write_statistics_lambda(self, value):
        self.stream.write(b"      <Statistics_lambda>%r</Statistics_lambda>\n" % value)

    def _write_statistics_entropy(self, value):
        self.stream.write(
            b"      <Statistics_entropy>%r</Statistics_entropy>\n" % value
        )

    def _start_iteration_hit(self):
        self.stream.write(b"<Hit>\n")

    def _end_iteration_hit(self):
        self.stream.write(b"</Hit>\n")

    def _write_hit_num(self, num):
        self.stream.write(b"  <Hit_num>%d</Hit_num>\n" % num)

    def _write_hit_id(self, hit_id):
        self.stream.write(b"  <Hit_id>%b</Hit_id>\n" % hit_id)

    def _write_hit_def(self, hit_def):
        self.stream.write(b"  <Hit_def>%b</Hit_def>\n" % hit_def)

    def _write_hit_accession(self, hit_accession):
        self.stream.write(b"  <Hit_accession>%b</Hit_accession>\n" % hit_accession)

    def _write_hit_len(self, hit_length):
        self.stream.write(b"  <Hit_len>%d</Hit_len>\n" % hit_length)

    def _start_hit_hsps(self):
        self.stream.write(b"  <Hit_hsps>\n")

    def _end_hit_hsps(self):
        self.stream.write(b"  </Hit_hsps>\n")

    def _start_hsp(self):
        self.stream.write(b"    <Hsp>\n")

    def _end_hsp(self):
        self.stream.write(b"    </Hsp>\n")

    def _write_hsp_num(self, num):
        self.stream.write(b"    <Hsp_num>%d</Hsp_num>\n" % num)

    def _write_hsp_bit_score(self, bit_score):
        self.stream.write(b"    <Hsp_bit-score>%b</Hsp_bit-score>\n" % bit_score)

    def _write_hsp_score(self, score):
        self.stream.write(b"    <Hsp_score>%d</Hsp_score>\n" % score)

    def _write_hsp_evalue(self, evalue):
        self.stream.write(b"    <Hsp_evalue>%b</Hsp_evalue>\n" % evalue)

    def _write_hsp_query_from(self, query_from):
        self.stream.write(b"     <Hsp_query-from>%d</Hsp_query-from>\n" % query_from)

    def _write_hsp_query_to(self, query_to):
        self.stream.write(b"     <Hsp_query-to>%d</Hsp_query-to>\n" % query_to)

    def _write_hsp_hit_from(self, hit_from):
        self.stream.write(b"     <Hsp_hit-from>%d</Hsp_hit-from>\n" % hit_from)

    def _write_hsp_hit_to(self, hit_to):
        self.stream.write(b"     <Hsp_hit-to>%d</Hsp_hit-to>\n" % hit_to)

    def _write_hsp_pattern_from(self, pattern_from):
        self.stream.write(
            b"     <Hsp_pattern-from>%d</Hsp_pattern-from>\n" % pattern_from
        )

    def _write_hsp_pattern_to(self, pattern_to):
        self.stream.write(b"     <Hsp_pattern-to>%d</Hsp_pattern-to>\n" % pattern_to)

    def _write_hsp_query_frame(self, query_frame):
        self.stream.write(b"     <Hsp_query-frame>%d</Hsp_query-frame>\n" % query_frame)

    def _write_hsp_hit_frame(self, hit_frame):
        self.stream.write(b"     <Hsp_hit-frame>%d</Hsp_hit-frame>\n" % hit_frame)

    def _write_hsp_identity(self, identity):
        self.stream.write(b"     <Hsp_identity>%d</Hsp_identity>\n" % identity)

    def _write_hsp_positive(self, positive):
        self.stream.write(b"     <Hsp_positive>%d</Hsp_positive>\n" % positive)

    def _write_hsp_gaps(self, gaps):
        self.stream.write(b"     <Hsp_gaps>%d</Hsp_gaps>\n" % gaps)

    def _write_hsp_align_len(self, align_len):
        self.stream.write(b"     <Hsp_align-len>%d</Hsp_align-len>\n" % align_len)

    def _write_hsp_density(self, density):
        self.stream.write(b"     <Hsp_density>%d</Hsp_density>\n" % density)

    def _write_hsp_qseq(self, qseq):
        self.stream.write(b"      <Hsp_qseq>%b</Hsp_qseq>\n" % qseq)

    def _write_hsp_hseq(self, hseq):
        self.stream.write(b"      <Hsp_hseq>%b</Hsp_hseq>\n" % hseq)

    def _write_hsp_midline(self, midline):
        self.stream.write(b"      <Hsp_midline>%b</Hsp_midline>\n" % midline)


class XML2Writer(BaseXMLWriter):
    """XML2 Writer."""

    def _start_blastoutput(self):
        self.stream.write(
            b"""\
<BlastXML2
    xmlns="http://www.ncbi.nlm.nih.gov"
    xmlns:xs="http://www.w3.org/2001/XMLSchema-instance"
    xs:schemaLocation="http://www.ncbi.nlm.nih.gov http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd"
>
<BlastOutput2>
  <report>
    <Report>
"""
        )

    def _end_blastoutput(self):
        self.stream.write(
            b"""\
    </Report>
  </report>
</BlastOutput2>
</BlastXML2>
"""
        )

    def _write_program(self, program):
        self.stream.write(b"      <program>%b</program>\n" % program)

    def _write_version(self, version):
        self.stream.write(b"      <version>%b</version>\n" % version)

    def _write_reference(self, reference):
        self.stream.write(b"      <reference>%b</reference>\n" % reference)

    def _write_db(self, db):
        self.stream.write(
            b"""\
      <search-target>
        <Target>
          <db>%b</db>
        </Target>
      </search-target>
"""
            % db
        )

    def _start_param(self):
        self.stream.write(
            b"""\
      <params>
        <Parameters>
"""
        )

    def _end_param(self):
        self.stream.write(
            b"""\
        </Parameters>
      </params>
"""
        )

    def _start_iterations(self):
        self.stream.write(
            b"""\
      <results>
        <Results>
"""
        )

    def _end_iterations(self):
        self.stream.write(
            b"""\
        </Results>
      </results>
"""
        )

    def _start_iteration(self):
        self.stream.write(
            b"""\
      <search>
        <Search>
"""
        )

    def _end_iteration(self):
        self.stream.write(
            b"""\
        </Search>
      </search>
"""
        )

    def _write_iteration_query_id(self, query_id):
        self.stream.write(
            b"""\
              <query-id>%b</query-id>
"""
            % query_id
        )

    def _write_iteration_query_def(self, query_def):
        self.stream.write(
            b"""\
              <query-title>%b</query-title>
"""
            % query_def
        )

    def _write_iteration_query_len(self, query_len):
        self.stream.write(
            b"""\
              <query-len>%d</query-len>
"""
            % query_len
        )

    def _write_iteration_query_masking(self, location):
        self.stream.write(
            b"""\
              <query-masking>
                <Range>
                  <from>%d</from>
                  <to>%d</to>
                </Range>
              </query-masking>
"""
            % (location.start + 1, location.end)
        )

    def _start_iteration_hits(self):
        self.stream.write(
            b"""\
              <hits>
"""
        )

    def _end_iteration_hits(self):
        self.stream.write(
            b"""\
              </hits>
"""
        )

    def _start_iteration_stat(self):
        self.stream.write(
            b"""\
              <stat>
"""
        )

    def _end_iteration_stat(self):
        self.stream.write(
            b"""\
              </stat>
"""
        )

    def _write_parameters_matrix(self, value):
        self.stream.write(
            b"""\
          <matrix>%b</matrix>
"""
            % value
        )

    def _write_parameters_expect(self, value):
        self.stream.write(
            b"""\
          <expect>%g</expect>
"""
            % value
        )

    def _write_parameters_include(self, value):
        self.stream.write(
            b"""\
          <include>%g</include>
"""
            % value
        )

    def _write_parameters_sc_match(self, value):
        self.stream.write(
            b"""\
          <sc-match>%d</sc-match>
"""
            % value
        )

    def _write_parameters_sc_mismatch(self, value):
        self.stream.write(
            b"""\
          <sc-mismatch>%d</sc-mismatch>
"""
            % value
        )

    def _write_parameters_gap_open(self, value):
        self.stream.write(
            b"""\
          <gap-open>%d</gap-open>
"""
            % value
        )

    def _write_parameters_gap_extend(self, value):
        self.stream.write(
            b"""\
          <gap-extend>%d</gap-extend>
"""
            % value
        )

    def _write_parameters_filter(self, value):
        self.stream.write(
            b"""\
          <filter>%b</filter>
"""
            % value
        )

    def _write_parameters_pattern(self, value):
        self.stream.write(
            b"""\
          <pattern>%b</pattern>
"""
            % value
        )

    def _write_parameters_entrez_query(self, value):
        self.stream.write(
            b"""\
          <entrez-query>%b</entrez-query>
"""
            % value
        )

    def _write_parameters_cbs(self, value):
        self.stream.write(
            b"""\
          <cbs>%d</cbs>
"""
            % value
        )

    def _write_parameters_query_gencode(self, value):
        self.stream.write(
            b"""\
          <query-gencode>%d</query-gencode>
"""
            % value
        )

    def _write_parameters_db_gencode(self, value):
        self.stream.write(
            b"""\
          <db-gencode>%d</db-gencode>
"""
            % value
        )

    def _write_parameters_bl2seq_mode(self, value):
        self.stream.write(
            b"""\
          <bl2seq-mode>%b</bl2seq-mode>
"""
            % value
        )

    def _write_statistics_db_num(self, db_num):
        self.stream.write(
            b"""\
                  <db-num>%d</db-num>
"""
            % db_num
        )

    def _write_statistics_db_len(self, db_len):
        self.stream.write(
            b"""\
                  <db-len>%d</db-len>
"""
            % db_len
        )

    def _write_statistics_hsp_len(self, hsp_len):
        self.stream.write(
            b"""\
                  <hsp-len>%d</hsp-len>
"""
            % hsp_len
        )

    def _write_statistics_eff_space(self, value):
        self.stream.write(
            b"""\
                  <eff-space>%s</eff-space>
"""
            % value
        )

    def _write_statistics_kappa(self, value):
        self.stream.write(
            b"""\
                  <kappa>%r</kappa>
"""
            % value
        )

    def _write_statistics_lambda(self, value):
        self.stream.write(
            b"""\
                  <lambda>%r</lambda>
"""
            % value
        )

    def _write_statistics_entropy(self, value):
        self.stream.write(
            b"""\
                  <entropy>%r</entropy>
"""
            % value
        )

    def _start_iteration_hit(self):
        self.stream.write(
            b"""\
                <Hit>
"""
        )

    def _end_iteration_hit(self):
        self.stream.write(
            b"""\
                </Hit>
"""
        )

    def _write_hit_num(self, num):
        self.stream.write(
            b"""\
                  <num>%d</num>
"""
            % num
        )

    def _start_hit_targets(self):
        self.stream.write(
            b"""\
                  <description>
"""
        )

    def _end_hit_targets(self):
        self.stream.write(
            b"""\
                  </description>
"""
        )

    def _start_hitdescr(self):
        self.stream.write(
            b"""\
                    <HitDescr>
"""
        )

    def _end_hitdescr(self):
        self.stream.write(
            b"""\
                    </HitDescr>
"""
        )

    def _write_hit_id(self, hit_id):
        self.stream.write(
            b"""\
                      <id>%b</id>
"""
            % hit_id
        )

    def _write_hit_def(self, hit_def):
        self.stream.write(
            b"""\
                      <title>%b</title>
"""
            % hit_def
        )

    def _write_hit_accession(self, hit_accession):
        self.stream.write(
            b"""\
                      <accession>%b</accession>
"""
            % hit_accession
        )

    def _write_hit_len(self, hit_length):
        self.stream.write(
            b"""\
                  <len>%d</len>
"""
            % hit_length
        )

    def _write_hit_taxid(self, taxid):
        self.stream.write(
            b"""\
                      <taxid>%d</taxid>
"""
            % taxid
        )

    def _write_hit_sciname(self, sciname):
        self.stream.write(
            b"""\
                      <sciname>%b</sciname>
"""
            % sciname
        )

    def _start_hit_hsps(self):
        self.stream.write(
            b"""\
                  <hsps>
"""
        )

    def _end_hit_hsps(self):
        self.stream.write(
            b"""\
                  </hsps>
"""
        )

    def _start_hsp(self):
        self.stream.write(
            b"""\
                    <Hsp>
"""
        )

    def _end_hsp(self):
        self.stream.write(
            b"""\
                    </Hsp>
"""
        )

    def _write_hsp_num(self, num):
        self.stream.write(
            b"""\
                      <num>%d</num>
"""
            % num
        )

    def _write_hsp_bit_score(self, bit_score):
        self.stream.write(
            b"""\
                      <bit-score>%b</bit-score>
"""
            % bit_score
        )

    def _write_hsp_score(self, score):
        self.stream.write(
            b"""\
                      <score>%d</score>
"""
            % score
        )

    def _write_hsp_evalue(self, evalue):
        self.stream.write(
            b"""\
                      <evalue>%b</evalue>
"""
            % evalue
        )

    def _write_hsp_query_from(self, query_from):
        self.stream.write(
            b"""\
                      <query-from>%d</query-from>
"""
            % query_from
        )

    def _write_hsp_query_to(self, query_to):
        self.stream.write(
            b"""\
                      <query-to>%d</query-to>
"""
            % query_to
        )

    def _write_hsp_hit_from(self, hit_from):
        self.stream.write(
            b"""\
                      <hit-from>%d</hit-from>
"""
            % hit_from
        )

    def _write_hsp_hit_to(self, hit_to):
        self.stream.write(
            b"""\
                      <hit-to>%d</hit-to>
"""
            % hit_to
        )

    def _write_hsp_pattern_from(self, pattern_from):
        self.stream.write(
            b"""\
                      <pattern-from>%d</pattern-from>
"""
            % pattern_from
        )

    def _write_hsp_pattern_to(self, pattern_to):
        self.stream.write(
            b"""\
                      <pattern-to>%d</pattern-to>
"""
            % pattern_to
        )

    def _write_hsp_query_frame(self, query_frame):
        self.stream.write(
            b"""\
                      <query-frame>%d</query-frame>
"""
            % query_frame
        )

    def _write_hsp_hit_frame(self, hit_frame):
        self.stream.write(
            b"""\
                      <hit-frame>%d</hit-frame>
"""
            % hit_frame
        )

    def _write_hsp_identity(self, identity):
        self.stream.write(
            b"""\
                      <identity>%d</identity>
"""
            % identity
        )

    def _write_hsp_positive(self, positive):
        self.stream.write(
            b"""\
                      <positive>%d</positive>
"""
            % positive
        )

    def _write_hsp_gaps(self, gaps):
        self.stream.write(
            b"""\
                      <gaps>%d</gaps>
"""
            % gaps
        )

    def _write_hsp_align_len(self, align_len):
        self.stream.write(
            b"""\
                      <align-len>%d</align-len>
"""
            % align_len
        )

    def _write_hsp_density(self, density):
        self.stream.write(
            b"""\
                      <density>%d</density>
"""
            % density
        )

    def _write_hsp_qseq(self, qseq):
        self.stream.write(
            b"""\
                      <qseq>%b</qseq>
"""
            % qseq
        )

    def _write_hsp_hseq(self, hseq):
        self.stream.write(
            b"""\
                      <hseq>%b</hseq>
"""
            % hseq
        )

    def _write_hsp_midline(self, midline):
        self.stream.write(
            b"""\
                      <midline>%b</midline>
"""
            % midline
        )
