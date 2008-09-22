# This is a Python module.
"""This module is DEPRECATED.

Andrew Dalke is no longer maintaining Martel or Bio.Mindy, and these modules
and associate ones like Bio.Std are now deprecated.  They are no longer
used in any of the current Biopython parsers, and are likely to be removed
in a future release.
"""

import warnings
warnings.warn("Martel and those parts of Biopython depending on it" \
              +" directly (such as Bio.Mindy and Bio.Std) are now" \
              +" deprecated, and will be removed in a future release of"\
              +" Biopython.  If you want to continue to use this code,"\
              +" please get in contact with the Biopython developers via"\
              +" the mailing lists to avoid its permanent removal from"\
              +" Biopython.", \
              DeprecationWarning)
# Standard Bioformats definitions

import Martel
Group = Martel.Group

namespace = "bioformat"
NS = namespace + ":"
XMLNS = "http://biopython.org/bioformat"

def _set_if_given(attrs, field, d, valid = None, convert = None):
    value = attrs.get(field)
    if value is not None:
        if valid is not None:
            if value not in valid:
                raise TypeError("%s (%r) must be one of %s" % \
                                (field, value, valid))
        if convert is None:
            d[field] = value
        else:
            d[field] = convert(value)

def _complain_if_given(attrs, name):
    if attrs.has_key(name) and attrs[name] is not None:
        raise NotImplementedError("Don't yet handle %r" % (name,))

def _must_have(expr, f):
    tag = f.tag
    if tag not in expr.group_names():
        raise TypeError(
            "group %r not present in the expression but is required" % \
            (tag,))

def _must_have_set(expr, sets):
    names = expr.group_names()
    for set in sets:
        for f in set:
            tag = f.tag
            if tag not in names:
                break
        else:
            return
    if len(sets) == 1:
        raise TypeError("missing required tags (need %s) in expression" %
                        [f.tag for f in sets[0]])
    lines = ["missing required tags in expression; must have one set from:"]
    for set in sets:
        lines.append( str( [t.tag for f in set] ) )
    s = "\n".join(lines)
    raise TypeError(s)

def _must_not_have(expr, f):
    f.tag
    if tag in expr.group_names():
        raise TypeError(
            "group %r present in the expression but is not allowed" % \
            (tag,))


# pre- Python 2.2 functions didn't allow attributes
def _f():
    pass
try:
    _f.x = 1
    _use_hack = 0
except AttributeError:
    _use_hack = 1
del _f

def _check_name(f, text):
    if text == "record": # XXX FIXME
        return
    assert NS + f.func_name == text, (NS + ":" + f.func_name, text)

def _check_attrs(attrs, names):
    for name in attrs.keys():
        if name not in names:
            raise TypeError("attr %r is not allowed here (valid terms: %s)" % \
                            (name, names))
    d = attrs.copy()
    for name in names:
        if not d.has_key(name):
            d[name] = None
    return d

if not _use_hack:
    def _settag(f, tag):
        _check_name(f, tag)
        f.tag = tag
else:
    # Convert the functions into callable objects
    class StdTerm:
        def __init__(self, func):
            self._func = func
        def __call__(self, *args, **kwargs):
            return self._func( *args, **kwargs)

    def _settag(f, tag):
        _check_name(f, tag)
        x = globals()[f.func_name] = StdTerm(f)
        x.tag = tag

################ identifier, description, and cross-references
def record(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("format",))
    d = {"xmlns:bioformat": XMLNS}
    _set_if_given(attrs, "format", d)
    return Group("record", expr, d) # XXX FIXME
_settag(record, "record") # XXX AND FIXME


def dbid(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("type", "style", "dbname"))
    d = {}
    _set_if_given(attrs, "type", d, ("primary", "accession", "secondary"))
    _set_if_given(attrs, "dbname", d)
    return Group(NS + "dbid", expr, d)
_settag(dbid, NS + "dbid")

def description_block(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("join",))
    _must_have(expr, description)
    d = {}
    _set_if_given(attrs, "join", d, ("english", "concat", "space", "newline"))
    return Group(NS + "description_block", expr, d)
_settag(description_block, NS + "description_block")

def description(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group(NS + "description", expr)
_settag(description, NS + "description")

def description_line(expr, attrs = {}):
    return description_block(description(expr, attrs))

def fast_dbxref(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("style",))
    d = {}
    _set_if_given(attrs, "style", d, ("sp-general", "sp-prosite", "sp-embl"))
    return Group(NS + "fast_dbxref", expr, d)

def dbxref(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("style",))
    _must_have(expr, dbxref_dbid)
    d = {}
    _complain_if_given(attrs, "style")
    return Group(NS + "dbxref", expr, d)
_settag(dbxref, NS + "dbxref")

def dbxref_dbname(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("style",))
    d = {}
    _set_if_given(attrs, "style", d)
    return Group(NS + "dbxref_dbname", expr, d)
_settag(dbxref_dbname, NS + "dbxref_dbname")

def dbxref_dbid(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("dbname", "type", "style", "negate"))
    d = {}
    _set_if_given(attrs, "dbname", d)
    _set_if_given(attrs, "type", d, ("primary", "accession", "secondary"))
    _complain_if_given(attrs, "style")
    _set_if_given(attrs, "negate", d, (0, 1), str)
    
    return Group(NS + "dbxref_dbid", expr, d)
_settag(dbxref_dbid, NS + "dbxref_dbid")

def dbxref_negate(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group(NS + "dbxref_negate", expr)
_settag(dbxref_negate, NS + "dbxref_negate")

##################### sequences

def _check_gapchar(s):
    if not ( ord(" ") <= ord(s) <= 126 ):
        raise TypeError("%r not allowed as a gap character" % (s,))
    return s

# What about three letter codes?
def sequence_block(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("alphabet", "gapchar", "remove_spaces"))
    _must_have(expr, sequence)
    d = {}
    _set_if_given(attrs, "alphabet", d,
                  ("iupac-protein", "iupac-dna", "iupac-rna",
                   "iupac-ambiguous-protein",
                   "iupac-ambiguous-dna",
                   "iupac-ambiguous-rna",
                   "protein", "dna", "rna", "unknown"))
    _set_if_given(attrs, "gapchar", d, convert = _check_gapchar)
    _set_if_given(attrs, "remove_spaces", d, (0, 1), str)
    return Group(NS + "sequence_block", expr, d)
_settag(sequence_block, NS + "sequence_block")

def sequence(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group(NS + "sequence", expr)
_settag(sequence, NS + "sequence")

def alphabet(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("alphabet",))
    d = {}
    _set_if_given(attrs, "alphabet", d,
                  ("iupac-protein", "iupac-dna", "iupac-rna",
                   "iupac-ambiguous-protein",
                   "iupac-ambiguous-dna",
                   "iupac-ambiguous-rna",
                   "protein", "dna", "rna", "nucleotide", "unknown"))
    return Group(NS + "alphabet", expr, d)
_settag(alphabet, NS + "alphabet")

    

############################## features

# In PIR

# FEATURE
#    1-25                #domain signal sequence #status predicted #label SIG\
#    26-737              #product procollagen-lysine 5-dioxygenase 2 #status
#                        predicted #label MAT\
#    63,209,297,365,522,
#    725                 #binding_site carbohydrate (Asn) (covalent) #status
#                        predicted

# The whole thing is a 'feature_block'

# One 'feature' is
#    26-737              #product procollagen-lysine 5-dioxygenase 2 #status
#                        predicted #label MAT\

# One 'feature_name' is "binding_site".

# An example of the feature_location_block and feature_block, which I
# will abbreviate as 'flb' and 'fl', is:
# <flb>   <fl>63,209,297,365,522,</fl>
#    <fl>725</fl>                 #binding_site carbohydrate ...

# PIR doesn't have a 'feature_description'

# Let:
#   fq = feature_qualifier
#   fqb = feature_qualifier
#   fqn = feature_qualifier_name
#   fqd = feature_qualifier_description
# then the text
#   
#    26-737              #product procollagen-lysine 5-dioxygenase 2 #status
#                        predicted #label MAT\
# 
# can be represented as (the rather tedious)
# 
#    26-737              <fqb><fq>#<fqn>product</fqn> <fqd>procollagen-\
# lysine 5-dioxygenase 2</fqd></fq> #<fq><fqn>status</fqn>
#                        <fqd>predicted</fqd> #<fq><fqn>label\
# </fqn> <fqd>MAT</fqd></fq>\</fqb>
#

# 'style' determines the namespace for the feature name
def feature_block(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("style", "location-style"))
    d = {}
    _set_if_given(attrs, "style", d)
    _set_if_given(attrs, "location-style", d)
    _must_have(expr, feature)
    return Group(NS + "feature_block", expr, d)
_settag(feature_block, NS + "feature_block")

def feature(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("location-style",))
    d = {}
    _set_if_given(attrs, "location-style", d)
    _must_have(expr, feature_name)
    _must_have_set(expr, [[feature_location],
                          [feature_location_start, feature_location_end]])
    return Group(NS + "feature", expr, d)
_settag(feature, NS + "feature")

def feature_name(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group(NS + "feature_name", expr)
_settag(feature_name, NS + "feature_name")

def feature_location(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group(NS + "feature_location", expr)
_settag(feature_location, NS + "feature_location")

def feature_location_start(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group(NS + "feature_location_start", expr)
_settag(feature_location_start, NS + "feature_location_start")

def feature_location_end(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group(NS + "feature_location_end", expr)
_settag(feature_location_end, NS + "feature_location_end")

def feature_description(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group(NS + "feature_description", expr)
_settag(feature_description, NS + "feature_description")


##def feature_qualifier_block(expr, attrs = {}):
##    attrs = _check_attrs(attrs, ())
##    _must_have(expr, feature_qualifier)
##    return Group(NS + "feature_qualifier_block", expr)
##_settag(feature_qualifier_block, NS + "feature_qualifier_block")

def feature_qualifier(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    _must_have(expr, feature_qualifier_name)
    return Group(NS + "feature_qualifier", expr)
_settag(feature_qualifier, NS + "feature_qualifier")

def feature_qualifier_name(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group(NS + "feature_qualifier_name", expr)
_settag(feature_qualifier_name, NS + "feature_qualifier_name")

def feature_qualifier_description(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group(NS + "feature_qualifier_description", expr)
_settag(feature_qualifier_description, NS + "feature_qualifier_description")


############ For homology searches

# "BLASTN", "BLASTP"
def application_name(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("app",))
    return Group("bioformat:application_name", expr, attrs)

# "2.0.11", "2.0a19MP-WashU"
def application_version(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:application_version", expr, attrs)

def search_header(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:search_header", expr, attrs)

def search_table(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:search_table", expr, attrs)

def search_table_description(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("bioformat:decode",))
    d = {"bioformat:decode": "strip"}
    _set_if_given(attrs, "bioformat:decode", d)
    return Group("bioformat:search_table_description", expr, d)

def search_table_value(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("name", "bioformat:decode"))
    return Group("bioformat:search_table_value", expr, attrs)

def search_table_entry(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:search_table_entry", expr, attrs)

def query_description_block(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("join-query",))
    d = {"join-query": "join|fixspaces"}
    _set_if_given(attrs, "join-query", d)
    return Group("bioformat:query_description_block", expr, d)

def query_description(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("bioformat:decode"))
    d = {}
    _set_if_given(attrs, "bioformat:decode", d)
    return Group("bioformat:query_description", expr, d)

def query_size(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:query_size", expr)

def database_name(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:database_name", expr, attrs)

def database_num_sequences(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("bioformat:decode",))
    return Group("bioformat:database_num_sequences", expr, attrs)

def database_num_letters(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("bioformat:decode",))
    return Group("bioformat:database_num_letters", expr, attrs)

def hit(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("join-description",))
    d = {"join-description": "join|fixspaces"}
    _set_if_given(attrs, "join-description", d)
    return Group("bioformat:hit", expr, d)

def hit_length(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:hit_length", expr, attrs)

def hit_description(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("bioformat:decode"))
    d = {}
    _set_if_given(attrs, "bioformat:decode", d)
    return Group("bioformat:hit_description", expr, d)

def hsp(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:hsp", expr, attrs)

def hsp_value(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("name", "bioformat:decode"))
    return Group("bioformat:hsp_value", expr, attrs)

def hsp_frame(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("which",))
    d = {}
    _set_if_given(attrs, "which", d, valid = ("query", "homology", "subject"))
    return Group("bioformat:hsp_frame", expr, d)

def hsp_strand(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("strand", "which"))
    d = {}
    _set_if_given(attrs, "which", d, valid = ("query", "homology", "subject"))
    _set_if_given(attrs, "strand", d, valid = ("+1", "0", "-1", ""))
    return Group("bioformat:hsp_strand", expr, d)

def hsp_seqalign_query_seq(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:hsp_seqalign_query_seq", expr, attrs)

def hsp_seqalign_homology_seq(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:hsp_seqalign_homology_seq", expr, attrs)

def hsp_seqalign_subject_seq(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:hsp_seqalign_subject_seq", expr, attrs)

def hsp_seqalign_query_leader(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:hsp_seqalign_query_leader", expr, attrs)
    

def hsp_seqalign_query_name(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:hsp_seqalign_query_name", expr, attrs)

def hsp_seqalign_subject_name(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:hsp_seqalign_subject_name", expr, attrs)

def hsp_seqalign(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:hsp_seqalign", expr, attrs)

def hsp_seqalign_query_start(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:hsp_seqalign_query_start", expr, attrs)

def hsp_seqalign_query_end(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:hsp_seqalign_query_end", expr, attrs)

def hsp_seqalign_subject_start(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:hsp_seqalign_subject_start", expr, attrs)

def hsp_seqalign_subject_end(expr, attrs = {}):
    attrs = _check_attrs(attrs, ())
    return Group("bioformat:hsp_seqalign_subject_end", expr, attrs)

def search_parameter(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("name", "bioformat:decode"))
    d = {}
    _set_if_given(attrs, "name", d)
    _set_if_given(attrs, "bioformat:decode", d)
    return Group("bioformat:search_parameter", expr, d)

def search_statistic(expr, attrs = {}):
    attrs = _check_attrs(attrs, ("name", "bioformat:decode"))
    d = {}
    _set_if_given(attrs, "name", d)
    _set_if_given(attrs, "bioformat:decode", d)
    return Group("bioformat:search_statistic", expr, d)

