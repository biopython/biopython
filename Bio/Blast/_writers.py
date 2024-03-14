import html
import codecs
from abc import ABC, abstractmethod


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

    def _write_hsp_query_frame(self, query_frame):
        self.stream.write(b"     <Hsp_query-frame>%d</Hsp_query-frame>\n" % query_frame)

    def _write_hsp_hit_frame(self, hit_frame):
        self.stream.write(b"     <Hsp_hit-frame>%d</Hsp_hit-frame>\n" % hit_frame)

    def _write_hsp_positive(self, positive):
        self.stream.write(b"     <Hsp_positive>%d</Hsp_positive>\n" % positive)

    def _write_hsp_identity(self, identity):
        self.stream.write(b"     <Hsp_identity>%d</Hsp_identity>\n" % identity)

    def _write_hsp_gaps(self, gaps):
        self.stream.write(b"     <Hsp_gaps>%d</Hsp_gaps>\n" % gaps)

    def _write_hsp_align_len(self, align_len):
        self.stream.write(b"     <Hsp_align-len>%d</Hsp_align-len>\n" % align_len)

    def _write_hsp_qseq(self, qseq):
        self.stream.write(b"      <Hsp_qseq>%b</Hsp_qseq>\n" % qseq)

    def _write_hsp_hseq(self, hseq):
        self.stream.write(b"      <Hsp_hseq>%b</Hsp_hseq>\n" % hseq)

    def _write_hsp_midline(self, midline):
        self.stream.write(b"      <Hsp_midline>%b</Hsp_midline>\n" % midline)

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
        positive = annotations["positive"]
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
        self._write_hsp_positive(positive)
        gaps = annotations.get("gaps")
        if gaps is not None:
            self._write_hsp_gaps(gaps)
        self._write_hsp_align_len(align_len)
        self._write_hsp_qseq(qseq.encode())
        self._write_hsp_hseq(hseq.encode())
        self._write_hsp_midline(midline.encode())
        self._end_hsp()

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

    def _write_hit(self, hit, query_length):
        stream = self.stream
        program = self._program
        target = hit.target
        target_length = len(target.seq)
        self._start_iteration_hit()
        self._write_hit_num(hit.num)
        self._write_hit_id(target.id.encode())
        self._write_hit_def(target.description.encode())
        self._write_hit_accession(target.name.encode())
        self._write_hit_len(target_length)
        self._start_hit_hsps()
        for hsp in hit:
            self._write_hsp(hsp, query_length, target_length)
        self._end_hit_hsps()
        self._end_iteration_hit()

    def _write_record(self, record):
        stream = self.stream
        self._start_iteration()
        self._write_iteration_num(record.num)
        query = record.query
        if query is None:
            query_length = None
        else:
            query_length = len(query.seq)
            self._write_iteration_query_id(query.id.encode())
            self._write_iteration_query_def(query.description.encode())
            self._write_iteration_query_len(query_length)
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

    def _write_xml_declaration(self):
        self.stream.write(b'<?xml version="1.0"?>\n')

    def _write_params(self, param):
        self._start_param()
        self.stream.write(b"    <Parameters>\n")
        value = param.get("matrix")
        if value is not None:
            self._write_parameters_matrix(value.encode())
        value = param.get("expect")
        if value is not None:
            self._write_parameters_expect(value)
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

    def _write_statistics(self, stat):
        self.stream.write(b"    <Statistics>\n")
        self.stream.write(
            b"      <Statistics_db-num>%d</Statistics_db-num>\n" % stat["db-num"]
        )
        self.stream.write(
            b"      <Statistics_db-len>%d</Statistics_db-len>\n" % stat["db-len"]
        )
        self.stream.write(
            b"      <Statistics_hsp-len>%d</Statistics_hsp-len>\n" % stat["hsp-len"]
        )
        self.stream.write(
            b"      <Statistics_eff-space>%s</Statistics_eff-space>\n"
            % str(stat["eff-space"]).encode()
        )
        self.stream.write(
            b"      <Statistics_kappa>%g</Statistics_kappa>\n" % stat["kappa"]
        )
        self.stream.write(
            b"      <Statistics_lambda>%g</Statistics_lambda>\n" % stat["lambda"]
        )
        self.stream.write(
            b"      <Statistics_entropy>%g</Statistics_entropy>\n" % stat["entropy"]
        )
        self.stream.write(b"    </Statistics>\n")

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
        self._write_query_id(records.query.id.encode())
        self._write_query_def(records.query.description.encode())
        self._write_query_len(len(records.query))
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

    @abstractmethod
    def _write_definition(self):
        return

    @abstractmethod
    def _start_blastoutput(self):
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
    def _write_query_id(self, query_id):
        return

    @abstractmethod
    def _write_query_def(self, query_def):
        return

    @abstractmethod
    def _write_query_len(self, query_len):
        return

    @abstractmethod
    def _start_param(self):
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

    @abstractmethod
    def _end_param(self):
        return

    @abstractmethod
    def _start_iterations(self):
        return

    @abstractmethod
    def _end_iterations(self):
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

    def _start_param(self):
        self.stream.write(b"  <BlastOutput_param>\n")

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

    def _start_iteration_stat(self):
        self.stream.write(b"  <Iteration_stat>\n")

    def _end_iteration_stat(self):
        self.stream.write(b"  </Iteration_stat>\n")
