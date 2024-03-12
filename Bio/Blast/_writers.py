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
        # fmt: off
        stream.write(f"""\
    <Hsp>
      <Hsp_num>{hsp.num}</Hsp_num>
      <Hsp_bit-score>{bit_score}</Hsp_bit-score>
      <Hsp_score>{hsp.score}</Hsp_score>
      <Hsp_evalue>{evalue}</Hsp_evalue>
      <Hsp_query-from>{query_from}</Hsp_query-from>
      <Hsp_query-to>{query_to}</Hsp_query-to>
      <Hsp_hit-from>{hit_from}</Hsp_hit-from>
      <Hsp_hit-to>{hit_to}</Hsp_hit-to>
      <Hsp_query-frame>{query_frame}</Hsp_query-frame>
      <Hsp_hit-frame>{hit_frame}</Hsp_hit-frame>
      <Hsp_identity>{identity}</Hsp_identity>
      <Hsp_positive>{positive}</Hsp_positive>
""".encode("UTF-8"))
        gaps = annotations.get("gaps")
        if gaps is not None:
            stream.write(b"""\
      <Hsp_gaps>%d</Hsp_gaps>
""" % gaps)
        stream.write(f"""\
      <Hsp_align-len>{align_len}</Hsp_align-len>
      <Hsp_qseq>{qseq}</Hsp_qseq>
      <Hsp_hseq>{hseq}</Hsp_hseq>
      <Hsp_midline>{midline}</Hsp_midline>
    </Hsp>
""".encode("UTF-8"))
        # fmt: on

    def _write_hit(self, hit, query_length):
        stream = self.stream
        program = self._program
        target = hit.target
        target_length = len(target.seq)
        # fmt: off
        stream.write(f"""\
<Hit>
  <Hit_num>{hit.num}</Hit_num>
  <Hit_id>{target.id}</Hit_id>
  <Hit_def>{target.description}</Hit_def>
  <Hit_accession>{target.name}</Hit_accession>
  <Hit_len>{target_length}</Hit_len>
  <Hit_hsps>
""".encode("UTF-8"))
        for hsp in hit:
            self._write_hsp(hsp, query_length, target_length)
        stream.write(b"""\
  </Hit_hsps>
</Hit>
""")
        # fmt: on

    def _write_record(self, record):
        stream = self.stream
        self._start_iteration()
        self._write_iteration_num(record.num)
        query = record.query
        if query is None:
            query_length = None
        else:
            query_length = len(query.seq)
            # fmt: off
            stream.write(f"""\
  <Iteration_query-ID>{query.id}</Iteration_query-ID>
  <Iteration_query-def>{query.description}</Iteration_query-def>
  <Iteration_query-len>{query_length}</Iteration_query-len>
""".encode("UTF-8"))
            # fmt: on
        # fmt: off
        stream.write(b"""\
<Iteration_hits>
""")
        # fmt: on
        for hit in record:
            self._write_hit(hit, query_length)
        # fmt: off
        stream.write(b"""\
</Iteration_hits>
""")
        # fmt: on
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
