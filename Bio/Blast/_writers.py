import html
import codecs


def _html_entity_replace(error):
    start = error.start
    end = error.end
    character = error.object[start:end]
    codepoint = ord(character)
    entity = f"&amp;{html.entities.codepoint2name[codepoint]};"
    return entity, end


codecs.register_error("htmlentityreplace", _html_entity_replace)


class XMLWriter:
    """XML Writer."""

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
        # fmt: off
        stream.write(f"""\
<Iteration>
  <Iteration_iter-num>{record.num}</Iteration_iter-num>
""".encode("UTF-8"))
        # fmt: on
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
            # fmt: off
            stream.write(f"""\
  <Iteration_stat>
    <Statistics>
      <Statistics_db-num>{stat["db-num"]}</Statistics_db-num>
      <Statistics_db-len>{stat["db-len"]}</Statistics_db-len>
      <Statistics_hsp-len>{stat["hsp-len"]}</Statistics_hsp-len>
      <Statistics_eff-space>{stat["eff-space"]}</Statistics_eff-space>
      <Statistics_kappa>{stat["kappa"]}</Statistics_kappa>
      <Statistics_lambda>{stat["lambda"]}</Statistics_lambda>
      <Statistics_entropy>{stat["entropy"]}</Statistics_entropy>
    </Statistics>
  </Iteration_stat>
""".encode("UTF-8"))
            # fmt: on
        # fmt: off
        stream.write(b"""\
</Iteration>
""")
        # fmt: on

    def write(self, records):
        """Write the records."""
        stream = self.stream
        count = 0
        program = records.program
        self._program = program
        # fmt: off
        stream.write(b"""\
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
"""
        )
        reference = html.escape(records.reference).encode("ASCII", "htmlentityreplace")
        reference = reference.decode()
        block = f"""\
  <BlastOutput_program>{program}</BlastOutput_program>
  <BlastOutput_version>{records.version}</BlastOutput_version>
  <BlastOutput_reference>{reference}</BlastOutput_reference>
  <BlastOutput_db>{records.db}</BlastOutput_db>
  <BlastOutput_query-ID>{records.query.id}</BlastOutput_query-ID>
  <BlastOutput_query-def>{records.query.description}</BlastOutput_query-def>
  <BlastOutput_query-len>{len(records.query)}</BlastOutput_query-len>
""".encode("UTF-8")
        stream.write(block)
        stream.write(b"""\
  <BlastOutput_param>
    <Parameters>
"""
        )
        param = records.param
        value = param.get("matrix")
        if value is not None:
            stream.write(b"""\
      <Parameters_matrix>%s</Parameters_matrix>
""" % value.encode("UTF-8"))
        value = param.get("expect")
        if value is not None:
            stream.write(b"""\
      <Parameters_expect>%g</Parameters_expect>
""" % value)
        value = param.get("include")
        if value is not None:
            stream.write(b"""\
      <Parameters_include>%g</Parameters_include>
""" % value)
        value = param.get("sc-match")
        if value is not None:
            stream.write(b"""\
      <Parameters_sc-match>%d</Parameters_sc-match>
""" % value)
        value = param.get("sc-mismatch")
        if value is not None:
            stream.write(b"""\
      <Parameters_sc-mismatch>%d</Parameters_sc-mismatch>
""" % value)
        value = param.get("gap-open")
        if value is not None:
            stream.write(b"""\
      <Parameters_gap-open>%d</Parameters_gap-open>
""" % value)
        value = param.get("gap-extend")
        if value is not None:
            stream.write(b"""\
      <Parameters_gap-extend>%d</Parameters_gap-extend>
""" % value)
        value = param.get("filter")
        if value is not None:
            stream.write(b"""\
      <Parameters_filter>%s</Parameters_filter>
""" % value.encode("UTF-8"))
        value = param.get("pattern")
        if value is not None:
            stream.write(b"""\
      <Parameters_pattern>%s</Parameters_pattern>
""" % value.encode("UTF-8"))
        value = param.get("entrez-query")
        if value is not None:
            stream.write(b"""\
      <Parameters_entrez-query>{value}</Parameters_entrez-query>
""" % value.encode("UTF-8"))
        stream.write(b"""\
    </Parameters>
  </BlastOutput_param>
<BlastOutput_iterations>
""")
        for record in records:
            self._write_record(record)
            count += 1
        stream.write(b"""\
</BlastOutput_iterations>
""")
        try:
            mbstat = records.mbstat
        except AttributeError:
            pass
        else:
            stream.write(f"""\
  <BlastOutput_mbstat>
    <Statistics>
      <Statistics_db-num>{mbstat["db-num"]}</Statistics_db-num>
      <Statistics_db-len>{mbstat["db-len"]}</Statistics_db-len>
      <Statistics_hsp-len>{mbstat["hsp-len"]}</Statistics_hsp-len>
      <Statistics_eff-space>{mbstat["eff-space"]}</Statistics_eff-space>
      <Statistics_kappa>{mbstat["kappa"]}</Statistics_kappa>
      <Statistics_lambda>{mbstat["lambda"]}</Statistics_lambda>
      <Statistics_entropy>{mbstat["entropy"]}</Statistics_entropy>
    </Statistics>
  </BlastOutput_mbstat>
""".encode("UTF-8"))
        stream.write(b"""\
</BlastOutput>
""")
        # fmt: on
        return count
