# Standard Content and Dispatch handlers for the Bioformat IO system

from xml.sax import handler
from Martel import Parser
from Bio import Std, Dispatch

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


class RecognizeHandler(handler.ContentHandler, handler.ErrorHandler):
    def __init__(self):
        self.recognized = 1
 
    def fatalError(self, exc):
        if isinstance(exc, Parser.ParserIncompleteException):
            pass
        else:
            self.recognized = 0
        raise exc
 
    error = fatalError

    def endElement(self, tag):
        if tag == "record":
            raise Parser.ParserException("we finished a record!")


def join_english(fields):
    if not fields:
        return ""
    s = fields[0]
    for field in fields[1:]:
        if s[-1:] == "-" and s[-3:-2] == "-":
            s = s + field
            continue
        if s.find(" ") == -1 and field.find(" ") == -1:
            s = s + field
            continue
        s = s + " " + field
    return s

def join_concat(fields):
    return fields.join("")

def join_space(fields):
    return fields.join(" ")

def join_newline(fields):
    return field.join("\n")
        

join_table = {
    "english": join_english,
    "concat": join_concat,
    "space": join_space,
    "newline": join_newline,
    }


class Handle_dbid(Dispatch.Callback):
    def start_dbid(self, tag, attrs):
        self.attrs = attrs
        self.save_characters()

    def end_dbid(self, tag):
        text = self.get_characters()
        self.callback(text, self.attrs)


class Handle_description(Dispatch.Callback):
    def start_description_block(self, tag, attrs):
        join = attrs.get("join", "english")
        self.join_fctn = join_table[join]
        self.descriptions = []
    def start_description(self, tag, attrs):
        self.save_characters()
    def end_description(self, tag):
        x = self.get_characters()
        self.descriptions.append(x)
    def end_description_block(self, tag):
        text = self.join_fctn(self.descriptions)
        self.descriptions = None
        self.callback(text)
    
class Handle_dbxref(Dispatch.Callback):
    def __init__(self, callback):
        Dispatch.Callback.__init__(self, callback)
        Dispatch.acquire_text(self, "dbxref_dbname", "dbname")
        
    def start_dbxref(self, tag, attrs):
        self.style = attrs.get("style", None)
        self.dbname = None
        self.ids = []
        self.negate = 0

    def start_dbxref_dbid(self, tag, attrs):
        dbname = attrs.get("dbname", None)
        type = attrs.get("type", None)
        self.ids.append( [None, dbname, type] )
        self.save_characters()

    def end_dbxref_dbid(self, tag):
        self.ids[-1][0] = self.get_characters()

    def end_dbxref(self, tag):
        negate = self.negate
        for id, dbname, type in self.ids:
            if dbname:
                style = None
            else:
                style = self.style
                dbname = self.dbname
            if dbname is None:
                raise TypeError("dbname not defined in the dbxref")
            self.callback( (dbname, id, type, negate) )

    def start_dbxref_negate(self, name, tag):
        self.negate = 1

class Handle_sequence(Dispatch.Callback):
    def __init__(self, callback):
        Dispatch.Callback.__init__(self, callback)
        Dispatch.acquire_append_text(self, "sequence", "sequence")

    def start_(self, tag, attrs):
        self.global_alphabet = None
        
    def start_sequence_block(self, tag, attrs):
        self.local_alphabet = attrs.get("alphabet", None)
        self.gapchar = attrs.get("gapchar", None)
        self.remove_spaces = attrs.get("remove_spaces", "1")
        self.sequence = []
        
    def end_sequence_block(self, tag):
        seq = "".join(self.sequence)
        if self.remove_spaces == "1":
            seq = seq.replace(" ", "")
        alphabet = self.local_alphabet or self.global_alphabet or "unknown"
        self.callback(alphabet, seq)

    def start_alphabet(self, tag, attrs):
        self.global_alphabet = attrs["alphabet"]

location_text_table = {
    "genbank": "".join,
    }

def convert_swissprot(start, end):
    return "Location(%r, %r)" % (start.strip(), end.strip())

location_pos_table = {
    None: lambda x, y: (int(x), int(y)),
    "swissprot": convert_swissprot,
    }

class Feature:
    def __init__(self, category, description, location, qualifiers):
        self.category = category
        self.description = description
        self.location = location
        self.qualifiers = qualifiers
    def __str__(self):
        return "Feature %r %r %s num_qualifiers = %d" % \
               (self.category, self.description, self.location,
                len(self.qualifiers))


class Handle_location(Dispatch.Callback):
    def __init__(self, callback, settings = {}):
        Dispatch.Callback.__init__(self, callback)
        self.settings = settings
        
        Dispatch.acquire_append_text(self, "_location", "location_text")
        Dispatch.acquire_text(self, "_location_start", "location_start")
        Dispatch.acquire_text(self, "_location_end", "location_end")

    # 'start_' and 'end_' are called by the 'feature' events
    def start_(self, tag, attrs):
        self.location_style = self.settings["location-style"]
        self.location_start = None
        self.location_end = None
        self.location_text = []

    def end_(self, tag):
        style = self.location_style
        if self.location_text:
            if self.location_start or self.location_end:
                raise TypeError("Cannot have both location text and start/end")
            function = location_text_table[style]
            location = function(self.location_text)
        else:
            function = location_pos_table[style]
            location = function(self.location_start, self.location_end)
        self.callback(location)
        

class Handle_feature_qualifier(Dispatch.Callback):
    def __init__(self, callback, settings):
        self.settings = settings
        Dispatch.Callback.__init__(self, callback)

        Dispatch.acquire_text(self, "feature_qualifier_name", "name")
        Dispatch.acquire_append_text(self, "feature_qualifier_description",
                                     "description")

    def start_feature_qualifier(self, tag, attrs):
        self.name = None
        self.description = []
        self.join = self.settings["qualifier-join"]

    def end_feature_qualifier(self, tag):
        self.callback( (self.name,
                        join_table[self.join](self.description)) )


class Handle_features(Dispatch.Callback):
    def __init__(self, callback):
        Dispatch.Callback.__init__(self, callback)
        self.settings = {}

        Dispatch.acquire_append_text(self, "feature_description",
                                     "description")

        Dispatch.acquire_text(self, tag = "feature_category",
                              attribute = "category")

        self.acquire(Handle_location(self.add_location, self.settings),
                     prefix = "feature")

        self.acquire(Handle_feature_qualifier(self.add_feature_qualifier,
                                              self.settings))

    def start_feature_block(self, tag, attrs):
        self.settings["join"] = attrs.get("join") or "english"
        self.settings["location-style"] = attrs.get("location-style") or \
                                          attrs.get("style") or \
                                          "swissprot"
        self.settings["qualifier-join"] = attrs.get("qualifier-join") or \
                                          attrs.get("join") or \
                                          "english"
        self.features = []

    def end_feature_block(self, tag):
        self.callback(self.features)
        self.features = None
        
    def start_feature(self, tag, attrs):
        self.category = None
        self.description = []
        self.location = None
        self.qualifiers = []

    def end_feature(self, tag):
        join = join_table[self.settings["join"]]
        self.features.append(Feature(self.category,
                                     join(self.description),
                                     self.location,
                                     self.qualifiers))
        
    def add_feature_qualifier(self, fq):
        self.qualifiers.append(fq)

    def add_location(self, location):
        self.location = location
