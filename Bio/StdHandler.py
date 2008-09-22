# Standard Content and Dispatch handlers for the Bioformat IO system
# This is a Python module.
"""This module is DEPRECATED.

Andrew Dalke is no longer maintaining Martel or Bio.Mindy, and these modules
and associate ones like Bio.StdHandler are now deprecated.  They are no longer
used in any of the current Biopython parsers, and are likely to be removed
in a future release.
"""

import warnings
warnings.warn("Martel and those parts of Biopython depending on it" \
              +" directly (such as Bio.Mindy and Bio.StdHandler) are now" \
              +" deprecated, and will be removed in a future release of"\
              +" Biopython.  If you want to continue to use this code,"\
              +" please get in contact with the Biopython developers via"\
              +" the mailing lists to avoid its permanent removal from"\
              +" Biopython.", \
              DeprecationWarning)

from xml.sax import handler
from Martel import Parser, Dispatch
from Bio import Std, Decode

###################################

# Helper functions to make functions

def add_int_handler(klass, tag, attrname):
    assert not hasattr(klass, "start_" +tag), "existing method exists"
    assert not hasattr(klass, "end_" +tag), "existing method exists"
    s = """if 1:
    def start(self, tag, attrs):
        self.save_characters()
    def end(self, tag):
        self.%s = int(self.get_characters())
""" % attrname
    d = {}
    exec s in d
    setattr(klass, "start_" + tag, d["start"])
    setattr(klass, "end_" + tag, d["end"])

def add_text_handler(klass, tag, attrname):
    assert not hasattr(klass, "start_" +tag), "existing method exists"
    assert not hasattr(klass, "end_" +tag), "existing method exists"
    s = """if 1:
    def start(self, tag, attrs):
        self.save_characters()
    def end(self, tag):
        self.%s = self.get_characters()
""" % attrname
    d = {}
    exec s in d
    setattr(klass, "start_" + tag, d["start"])
    setattr(klass, "end_" + tag, d["end"])

def add_text_dict_handler(klass, tag, attrname, key):
    assert not hasattr(klass, "start_" +tag), "existing method exists"
    assert not hasattr(klass, "end_" +tag), "existing method exists"
    s = """if 1:
    def start(self, tag, attrs):
        self.save_characters()
    def end(self, tag):
        self.%s["%s"] = self.get_characters()
""" % (attrname, key)
    d = {}
    exec s in d
    setattr(klass, "start_" + tag, d["start"])
    setattr(klass, "end_" + tag, d["end"])

def add_text_decode_handler(klass, tag, attrname):
    assert not hasattr(klass, "start_" +tag), "existing method exists"
    assert not hasattr(klass, "end_" +tag), "existing method exists"
    s = """if 1:
    def start(self, tag, attrs):
        self.save_characters()
        self._decode_%s = attrs.get("bioformat:decode", None)
    def end(self, tag):
        if self._decode_%s is not None:
            s = Decode.make_decoder(self._decode_%s)(s)
        self.%s = self.get_characters()
""" % (tag, tag, tag, attrname)
    d = {"Decode": Decode}
    exec s in d
    setattr(klass, "start_" + tag, d["start"])
    setattr(klass, "end_" + tag, d["end"])

def add_first_text_handler(klass, tag, attrname):
    assert not hasattr(klass, "start_" +tag), "existing method exists"
    assert not hasattr(klass, "end_" +tag), "existing method exists"
    s = """if 1:
    def start(self, tag, attrs):
        if self.%s is None:
            self.save_characters()
    def end(self, tag):
        if self.%s is None:
            self.%s = self.get_characters()
""" % (attrname, attrname, attrname)
    d = {}
    exec s in d
    setattr(klass, "start_" + tag, d["start"])
    setattr(klass, "end_" + tag, d["end"])

def add_text_block_handler(klass, tag, joinattr, defaultjoin, attrname):
    assert not hasattr(klass, "start_" + tag), "existing method exists"
    assert not hasattr(klass, "end_" + tag), "existing method exists"
    assert not hasattr(klass, "start_"+tag+"_block"), "existing method exists"
    assert not hasattr(klass, "end_" +tag+"_block"), "existing method exists"
    s = """if 1:
    def start_block(self, tag, attrs):
        self._%(tag)s_join_func = Decode.make_decoder(attrs.get(%(joinattr)r, %(defaultjoin)r))
        self._%(tag)s_lines = []
    def end_block(self, tag):
        self.%(attrname)s = self._%(tag)s_join_func(self._%(tag)s_lines)
    def start(self, tag, attrs):
        self.save_characters()
    def end(self, tag):
        self._%(tag)s_lines.append(self.get_characters())
""" % locals()
    d = {"Decode": Decode}
    exec s in d
    setattr(klass, "start_" + tag, d["start"])
    setattr(klass, "end_" + tag, d["end"])
    setattr(klass, "start_" + tag + "_block", d["start_block"])
    setattr(klass, "end_" + tag + "_block", d["end_block"])

def add_value_handler(klass, tag, attrname):
    assert not hasattr(klass, "start_" +tag), "existing method exists"
    assert not hasattr(klass, "end_" +tag), "existing method exists"
    s = """if 1:
    def start(self, tag, attrs):
        self._%(tag)s_name = attrs["name"]
        self._%(tag)s_decode = attrs.get("bioformat:decode", None)
        self.save_characters()
    def end(self, tag):
        s = self.get_characters()
        if self._%(tag)s_decode is not None:
            s = Decode.make_decoder(self._%(tag)s_decode)(s)
        self.%(attrname)s[self._%(tag)s_name] = s
""" % locals()
    d = {"Decode": Decode}
    exec s in d
    setattr(klass, "start_" + tag, d["start"])
    setattr(klass, "end_" + tag, d["end"])

    
#################################

class ConvertHandler(handler.ContentHandler):
    """Used to read records and produce output"""
    def __init__(self, record_builder, writer, record_tag = "record"):
        handler.ContentHandler.__init__(self)
        self.record_builder = record_builder
        self.writer = writer
        self.record_tag = record_tag

    def startDocument(self):
        self.inside_record = 0
        self.characters = self.ignore_characters
        
    def startElement(self, tag, attrs):
        if self.inside_record:
            self.record_builder.startElement(tag, attrs)
        elif tag == self.record_tag:
            self.record_builder.startDocument()
            self.inside_record = 1
            self.characters = self.record_builder.characters
            self.record_builder.startElement(tag, attrs)

    def endElement(self, tag):
        if self.inside_record:
            self.record_builder.endElement(tag)
            if tag == self.record_tag:
                self.record_builder.endDocument()
                self.writer.write(self.record_builder.document)
                self.inside_record = 0
                self.characters = self.ignore_characters

    def ignore_characters(self, s):
        pass

class ConvertDispatchHandler(Dispatch.Dispatcher):
    """Used to read records and produce output through a Dispatcher"""
    def __init__(self, record_builder, writer, record_tag = "record"):
        setattr(self, "end_" + record_tag, self.write_record)
        Dispatch.Dispatcher.__init__(self,
                                     remap = {record_tag: "bioformat:"}
                                     )
        self.acquire(record_builder)
        self.record_builder = record_builder
        self.writer = writer
        self.record_tag = record_tag
    def write_record(self, tag):
        self.writer.write(self.record_builder.document)



class RecognizeHandler(handler.ContentHandler, handler.ErrorHandler):
    def __init__(self):
        self.recognized = 1
        self.exc = None
 
    def fatalError(self, exc):
        if isinstance(exc, Parser.ParserIncompleteException):
            pass
        else:
            self.recognized = 0
            self.exc = exc
        raise exc
 
    error = fatalError

    def endElement(self, tag):
        if tag == "record":
            raise Parser.ParserException("we finished a record!")



class Handle_dbid(Dispatch.Callback):
    def start_dbid(self, tag, attrs):
        self.attrs = attrs
        self.save_characters()

    def end_dbid(self, tag):
        text = self.get_characters()
        self.callback(text, self.attrs)


class Handle_description(Dispatch.Callback):
    def start_description_block(self, tag, attrs):
        j = attrs.get("join", None)
        if j is None:
            self.join_fctn = Decode.join_fixspaces
        else:
            self.join_fctn = Decode.make_typechecked_decoder(j, list, str)
        self.descriptions = []
    def start_description(self, tag, attrs):
        self.save_characters()
    def end_description(self, tag):
        x = self.get_characters()
        self.descriptions.append(x)
    def end_description_block(self, tag):
        self.callback(self.join_fctn(self.descriptions))

#### There can be multiple dbxref_dbids in a dbxref
# DR   EMBL; X64411; CAA45756.1; -.
#    <dbxref><..dbname style="swiss">EMBL</..dbname>
#                        <dbid type="primary">X64411</dbid>
#                        <dbid type="accession">CAA45756.1</dbid>
#    </dbxref>
###
# DR   P35156, YPUI_BACSU, F;
#   <dbxref><dbid type="primary" dbname="sprot">P35156</dbid>
#           <dbid type="accession" dbname="sprot">YPUI_BACSU</dbid>
#           <negate/>
#    </dbxref>

def _fixup_sp_pattern(exp):
    import re
    import Martel
    exp = Martel.select_names(exp, (Std.dbxref_dbname.tag,Std.dbxref_dbid.tag))
                               
    e = exp._find_groups(Std.dbxref_dbname.tag)
    assert len(e) == 1
    e = e[0]
    e.name = "dbname"
    dbstyle = e.attrs["style"]
    e.attrs = {}
    e = exp._find_groups(Std.dbxref_dbid.tag)
    assert len(e) == 2
    e[0].name = "primary_dbid"
    primary_type = e[0].attrs["type"]
    e[0].attrs = {}
    e[1].name = "secondary_dbid"
    secondary_type = e[1].attrs["type"]
    e[1].attrs = {}
    pattern = str(exp) + "$"
    pat = re.compile(pattern)
    return pat, dbstyle, primary_type, secondary_type

# Turns out these 'fast' versions speed up the dbxref code by about
# a factor of 2.

# DR   PIR; S08427; S08427.
_fast_dbxref_sp_general_data = None
def _fast_dbxref_sp_general(s):
    global _fast_dbxref_sp_general_data
    if _fast_dbxref_sp_general_data is None:
        from Bio.expressions.swissprot import sprot38
        _fast_dbxref_sp_general_data = _fixup_sp_pattern(
                                                    sprot38.real_DR_general)

    pat, dbstyle, primary_type, secondary_type = _fast_dbxref_sp_general_data

    m = pat.match(s)
    assert m is not None, "Ill-formated sp-general dxbref: %r" % s
    return (
        (dbstyle, m.group("dbname"), primary_type,
                                     m.group("primary_dbid"), 0),
        (dbstyle, m.group("dbname"), secondary_type,
                                     m.group("secondary_dbid"), 0)
        )

# DR   PFAM; PF01018; GTP1_OBG; 1.
# DR   PROSITE; PS00905; GTP1_OBG; 1.

_fast_dbxref_sp_prosite_data = None
def _fast_dbxref_sp_prosite(s):
    global _fast_dbxref_sp_prosite_data

    if _fast_dbxref_sp_prosite_data is None:
        from Bio.expressions.swissprot import sprot38
        _fast_dbxref_sp_prosite_data = _fixup_sp_pattern(
                                                    sprot38.real_DR_prosite)

    pat, dbstyle, primary_type, secondary_type = _fast_dbxref_sp_prosite_data
    m = pat.match(s)
    assert m is not None, "Ill-formated sp-prosite dxbref: %r" % s
    return (
        (dbstyle, m.group("dbname"), primary_type,
                                     m.group("primary_dbid"), 0),
        (dbstyle, m.group("dbname"), secondary_type,
                                     m.group("secondary_dbid"), 0)
        )
    

# DR   EMBL; M36407; AAA33110.1; -.
_fast_dbxref_sp_embl_data = None
def _fast_dbxref_sp_embl(s):
    global _fast_dbxref_sp_embl_data

    if _fast_dbxref_sp_embl_data is None:
        from Bio.expressions.swissprot import sprot38
        _fast_dbxref_sp_embl_data = _fixup_sp_pattern(
                                                    sprot38.real_DR_embl)

    pat, dbstyle, primary_type, secondary_type = _fast_dbxref_sp_embl_data
    m = pat.match(s)
    assert m is not None, "Ill-formated sp-embl dxbref: %r" % s
    return (
        (dbstyle, m.group("dbname"), primary_type,
                                     m.group("primary_dbid"), 0),
        (dbstyle, m.group("dbname"), secondary_type,
                                     m.group("secondary_dbid"), 0)
        )

_fast_dbxref_parser_table = {
    "sp-general": _fast_dbxref_sp_general,
    "sp-prosite": _fast_dbxref_sp_prosite,
    "sp-embl": _fast_dbxref_sp_embl,
}

class Handle_dbxref(Dispatch.Callback):
    def __init__(self, callback):
        Dispatch.Callback.__init__(self, callback)
        self.supported_features.append("fast-sp-dbxref")
        self.slow_callback = self.callback
    def start_dbxref(self, tag, attrs):
        self.negate = 0
        self.dbname = None
        self.dbids = []
        self.info = []

    def start_dbxref_dbname(self, tag, attrs):
        assert self.dbname is None, "cannot set the dbname twice"
        self.dbname_style = attrs.get("style", "unknown")
        self.save_characters()
    def end_dbxref_dbname(self, tag):
        self.dbname = self.get_characters()

    def start_dbxref_dbid(self, tag, attrs):
        d = attrs.get("dbname", None)
        if d is None:
            assert self.dbname is not None, "must set the dbname"
            self.info.append( (self.dbname_style, self.dbname,
                               attrs.get("type", "primary")) )
        else:
            self.info.append( ("bioformat", d,
                               attrs.get("type", "primary")) )
        self.save_characters()

    def end_dbxref_dbid(self, tag):
        self.dbids.append( self.get_characters())

    def start_dbxref_negate(self, tag, attrs):
        self.negate = 1

    def end_dbxref(self, tag):
        cb = self.slow_callback
        if cb is None:
            return
        negate = self.negate
        for ( (dbname_style, dbname, idtype), dbid) in zip(self.info,
                                                           self.dbids):
            self.slow_callback(dbname_style, dbname, idtype, dbid, negate)

    def start_fast_dbxref(self, tag, attrs):
        style = attrs["style"]
        self._fast_parser = _fast_dbxref_parser_table[style]
        self.save_characters()
        self.slow_callback = None
    def end_fast_dbxref(self, tag):
        for info in self._fast_parser(self.get_characters()):
            self.callback(*info)
        self.slow_callback = self.callback

##################
class Handle_sequence(Dispatch.Callback):
    global_alphabet = None
    def start_(self, tag, attrs):
        self.global_alphabet = None
        
    def start_sequence_block(self, tag, attrs):
        self.local_alphabet = attrs.get("alphabet", None)
        self.gapchar = attrs.get("gapchar", None)
        self.stopchar = attrs.get("stopchar", None)
        j = attrs.get("join", None)
        if j is not None:
            self.join_func = Decode.make_typechecked_decoder(j, list, str)
        else:
            self.join_func = None
        self.sequences = []
        
    def end_sequence_block(self, tag):
        f = self.join_func
        if f is not None:
            seq = self.f(self.sequences)
        else:
            seq = "".join(self.sequences).replace(" ", "")
        alphabet = self.local_alphabet or self.global_alphabet or "unknown"
        self.callback( (alphabet, seq, self.gapchar, self.stopchar) )

    def start_alphabet(self, tag, attrs):
        self.global_alphabet = attrs["alphabet"]

    def start_sequence(self, tag, attrs):
        self.save_characters()
    def end_sequence(self, tag):
        self.sequences.append(self.get_characters())

class Feature:
    def __init__(self, name, description, location, qualifiers):
        self.name = name
        self.description = description
        self.location = location
        self.qualifiers = qualifiers
    def __str__(self):
        return "Feature %r %r %s num_qualifiers = %d" % \
               (self.name, self.description, self.location,
                len(self.qualifiers))


class Handle_feature_location(Dispatch.Callback):
    def __init__(self, callback, settings = {}):
        Dispatch.Callback.__init__(self, callback)
        self.settings = settings
        
    def start_feature(self, tag, attrs):
        self.location_style = attrs.get("location-style",
                                        self.settings["location-style"])
        j = attrs.get("join-feature", None)
        if j is None:
            self.text_join_func = "".join
        else:
            self.text_join_func = Decode.make_typechecked_decoder(j, list, str)

        self.location_start = None
        self.location_end = None
        self.text_lines = []

    def end_feature(self, tag):
        if self.location_start or self.location_end:
            if self.text_lines:
                raise TypeError("Cannot have both location text and start/end")
            self.callback(self.location_style,
                          (self.location_start, self.location_end))
        else:
            self.callback(self.location_style,
                          (self.text_join_func(self.text_lines), None))
    
    def start_feature_location(self, tag, attrs):
        self.save_characters()
    def end_feature_location(self, tag):
        self.text_lines.append(self.get_characters())
        
add_text_handler(Handle_feature_location, "feature_location_start",
                 "location_start")
add_text_handler(Handle_feature_location, "feature_location_end",
                 "location_end")

##################################

class Handle_feature_qualifier(Dispatch.Callback):
    def __init__(self, callback, settings):
        self.settings = settings
        Dispatch.Callback.__init__(self, callback)

    def start_feature_qualifier(self, tag, attrs):
        self.name = None
        self.description = []
        qj = attrs.get("join-qualifier", None)
        if qj is None:
            self.join = self.settings["qualifier_join_func"]
        else:
            self.join = Decode.make_typechecked_decoder(qj, list, str)

    def end_feature_qualifier(self, tag):
        self.callback(self.name, self.join(self.description))

    def start_feature_qualifier_description(self, tag, attrs):
        self.save_characters()
    def end_feature_qualifier_description(self, tag):
        self.description.append(self.get_characters())

add_text_handler(Handle_feature_qualifier, "feature_qualifier_name", "name")

####################

class Handle_features(Dispatch.Callback):
    def __init__(self, callback):
        Dispatch.Callback.__init__(self, callback)
        self.settings = {}

        self.acquire(Handle_feature_location(self.add_location, self.settings))

        self.acquire(Handle_feature_qualifier(self.add_feature_qualifier,
                                              self.settings))

    def start_feature_block(self, tag, attrs):
        jf = attrs.get("join-description", None)
        if jf is None:
            self.join_feature_description = Decode.join_fixspaces
        else:
            self.join_feature_description = Decode.make_typechecked_decoder(
                jf, list, str)

        self.settings["location-style"] = attrs.get("location-style", None)

        jq = attrs.get("join-qualifier", None)
        if jq is None:
            self.settings["qualifier_join_func"] = Decode.join_fixspaces
        else:
            self.settings["qualifier_join_func"] = \
                          Decode.make_typechecked_decoder(jq, list, str)
        self.features = []

    def end_feature_block(self, tag):
        self.callback(self.features)
        self.features = None
        
    def start_feature(self, tag, attrs):
        self.name = None
        self.description = []
        self.location = None
        self.qualifiers = []

    def start_feature_description(self, tag, attrs):
        self.save_characters()
    def end_feature_description(self, tag):
        self.description.append(self.get_characters())

    def end_feature(self, tag):
        self.features.append(Feature(
            self.name,
            self.join_feature_description(self.description),
            self.location,
            self.qualifiers))
        
    def add_feature_qualifier(self, name, description):
        self.qualifiers.append((name, description))

    def add_location(self, style, location_info):
        self.location = (style, location_info)

add_text_handler(Handle_features, "feature_name", "name")


############## Search handlers

class Handle_hsp_seqalign(Dispatch.Callback):
    def start_hsp(self, tag, attrs):
        self.query_name = None     # "Query"
        self.subject_name = None   # "Sbjct"
        
        self.query_seq = ""        # the actual text of the sequence
        self.homology_seq = ""
        self.subject_seq = ""
        
        self.query_start_loc = None
        self.query_end_loc = None
        
        self.subject_start_loc = None
        self.subject_end_loc = None

    def end_hsp(self, tag):
        self.callback(self)

    def start_hsp_seqalign(self, tag, attrs):
        self.sub_leader = None

    def start_hsp_seqalign_query_seq(self, tag, attrs):
        self.save_characters()
    def end_hsp_seqalign_query_seq(self, tag):
        s = self.get_characters()
        self.query_seq += s
        self.sub_query_seq_len = len(s)

    def start_hsp_seqalign_homology_seq(self, tag, attrs):
        self.save_characters()
    def end_hsp_seqalign_homology_seq(self, tag):
        query_leader = self.leader_size
        query_seq_len = self.sub_query_seq_len
        line = self.get_characters()
        s = line[query_leader:query_leader+query_seq_len]
        assert len(s) == query_seq_len, (len(s), query_seq_len, line)
        self.homology_seq += s

    def start_hsp_seqalign_subject_seq(self, tag, attrs):
        self.save_characters()
    def end_hsp_seqalign_subject_seq(self, tag):
        self.subject_seq += self.get_characters()
    
    def start_hsp_seqalign_query_leader(self, tag, attrs):
        self.save_characters()
    def end_hsp_seqalign_query_leader(self, tag):
        self.leader_size = len(self.get_characters())

add_first_text_handler(Handle_hsp_seqalign, "hsp_seqalign_query_name",
                         "query_name")

add_first_text_handler(Handle_hsp_seqalign, "hsp_seqalign_subject_name",
                         "subject_name")

add_first_text_handler(Handle_hsp_seqalign, "hsp_seqalign_query_start",
                         "query_start_loc")
add_text_handler(Handle_hsp_seqalign, "hsp_seqalign_query_end",
                 "query_end_loc")

add_first_text_handler(Handle_hsp_seqalign, "hsp_seqalign_subject_start",
                         "subject_start_loc")
add_text_handler(Handle_hsp_seqalign, "hsp_seqalign_subject_end",
                 "subject_end_loc")




#############################

class Handle_hsp(Dispatch.Callback):
    def __init__(self, callback):
        Dispatch.Callback.__init__(self, callback)
        self.acquire(Handle_hsp_seqalign(self.add_hsp_seqs))

    def start_hsp(self, tag, attrs):
        self.hsp_values = {}      # expect, p, identities, ...
        self.strands = {}
        self.frames = {}

    def end_hsp(self, tag):
        self.callback(self.hsp_values,
                      self.hsp_info,
                      self.strands, self.frames,
                      )
        
    def start_hsp_strand(self, tag, attrs):
        self.strands[attrs["which"]] = attrs["strand"]

    def start_hsp_frame(self, tag, attrs):
        self.getting_frame = attrs["which"]
        self.save_characters()

    def end_hsp_frame(self, tag):
        self.frames[self.getting_frame] = self.get_characters()
        self.getting_frame = None

    def add_hsp_seqs(self, hsp_info):
        self.hsp_info = hsp_info

    def start_hsp_value(self, tag, attrs):
        self.value_convert = attrs.get("bioformat:decode", None)
        self.value_name = attrs["name"]
        self.save_characters()

    def end_hsp_value(self, tag):
        s = self.get_characters()
        if self.value_name is not None:
            if self.value_name == "float":
                s = float(s)
            else:
                s = Decode.make_decoder(self.value_convert)(s)
        self.hsp_values[self.value_name] = s

#############################


class Handle_search_table(Dispatch.Callback):
    def start_search_table_value(self, tag, attrs):
        self.value_name = attrs["name"]
        self.value_decode = attrs.get("bioformat:decode", None)
        self.save_characters()
    def end_search_table_value(self, tag):
        s = self.get_characters()
        if self.value_decode is not None:
            x = self.value_decode
            if x == "int":
                s = int(s)
            elif x == "float":
                s = float(s)
            else:
                s = Decode.make_decoder(x)(s)
        self.values[self.value_name] = s

    def start_search_table(self, tag, attrs):
        self.data = []
    def end_search_table(self, tag):
        self.callback(self.data)
        self.data = None

    def start_search_table_entry(self, tag, attrs):
        self.description = None
        self.values = {}
        
    def end_search_table_entry(self, tag):
        self.data.append( (self.description, self.values) )
        self.description = self.values = None

add_text_handler(Handle_search_table, "search_table_description",
                 "description")

#############################
    
class Handle_search_header(Dispatch.Callback):
    def start_(self, tag, attrs):
        self.dict = {}
        self.query_description = None

    def end_search_header(self, tag):
        d = self.dict
        d["query_description"] = self.query_description
        self.callback(d)

add_text_block_handler(Handle_search_header, "query_description",
                       "join-query", "join|fixspaces", "query_description")

add_text_dict_handler(Handle_search_header, "application_name",
                      "dict", "appname")
add_text_dict_handler(Handle_search_header, "application_version",
                      "dict", "appversion")
add_text_dict_handler(Handle_search_header, "database_name",
                      "dict", "dbname")
add_text_dict_handler(Handle_search_header, "database_num_sequences",
                      "dict", "db_num_sequences")
add_text_dict_handler(Handle_search_header, "database_num_letters",
                      "dict", "db_num_letters")
add_text_dict_handler(Handle_search_header, "query_size",
                      "dict", "query_size")


#############################

class Handle_search_info(Dispatch.Callback):
    def start_(self, tag, attrs):
        self.parameters = {}
        self.statistics = {}
        
    def end_(self, tag):
        self.callback(self.parameters, self.statistics)

add_value_handler(Handle_search_info, "search_parameter", "parameters")
add_value_handler(Handle_search_info, "search_statistic", "statistics")
