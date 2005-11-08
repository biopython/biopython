import sys, urllib
from xml.sax import saxutils
import StringIO

from Bio.EUtils import ReseekFile
from Bio.config.FormatRegistry import FormatObject

class FormatIOIterator:
    def __init__(self, obj):
        self.obj = obj
        self._n = 0
    def next(self):
        x = self.obj.next()
        if x is not None:
            self._n += 1
            return x.document
        return x
    def __getitem__(self, i):
        assert i == self._n, "forward iteration only"
        x = self.next()
        if x is None:
            raise IndexError(i)
        return x

def _get_selected_names(builder, format):
    if not hasattr(builder, "uses_tags"):
        return None
    
    select_names = tuple(builder.uses_tags()) + ("record",)
    if not hasattr(builder, "get_supported_features"):
        return select_names
    
    d = {}
    for name in select_names:
        d[name] = 1
    supported_features = builder.get_supported_features()

    for name, remove_tags in format.expression.features():
        if name in supported_features:
            for tag in remove_tags:
                if d.has_key(tag):
                    del d[tag]
    return d.keys()

class FormatIO:
    def __init__(self, name,
                 default_input_format = None,
                 default_output_format = None,
                 abbrev = None,
                 registery = None):
        if abbrev is None:
            abbrev = name
        if registery is None:
            import Bio
            registery = Bio.formats
        
        self.name = name
        self.abbrev = abbrev
        self.default_input_format = default_input_format
        self.default_output_format = default_output_format
        self.registery = registery

    def _get_file_format(self, format, source):
        # By construction, this is the lowest we can go, so don't try
        if isinstance(format, FormatObject):
            return format, source

        source = saxutils.prepare_input_source(source)
        infile = source.getCharacterStream() or source.getByteStream()
        
        is_reseek = 0
        try:
            infile.tell()
        except (AttributeError, IOError):
            infile = ReseekFile.ReseekFile(infile)
            is_reseek = 1

        format = format.identifyFile(infile)

        if is_reseek:
            infile.nobuffer()

        # returned file could be a ReseekFile!
        return format, infile
        
    def readFile(self, infile, format = None, builder = None, debug_level = 0):
        if format is None:
            format = self.default_input_format
            if format is None:
                raise ValueError, "No format specified"
        format = self.registery.normalize(format)
        if builder is None:
            builder = self.registery.find_builder(format, self)
        format, infile = self._get_file_format(format, infile)
            
        if format is None:
            raise TypeError("Could not determine file type")

        select_names = _get_selected_names(builder, format)

        if format.multirecord == 1:
            iterator = format.make_iterator("record",
                                            select_names = select_names,
                                            debug_level = debug_level)
            return FormatIOIterator(iterator.iterateFile(infile, builder))
        elif format.multirecord == 0:
            parser = format.make_parser(select_names = select_names,
                                        debug_level = debug_level)
            parser.setContentHandler(builder)
            parser.parseFile(infile)
            return builder.document
        else:
            raise AssertionError(format.multirecord)

    def readString(self, s, *args, **keywds):
        infile = StringIO.StringIO(s)
        return self.readFile(infile, *args, **keywds)
                  
    def read(self, systemID, format = None, builder = None, debug_level = 0):
        if isinstance(systemID, type("")):
            # URL
            infile = ReseekFile.ReseekFile(urllib.urlopen(systemID))
        else:
            infile = systemID
        return self.readFile(infile, format, builder, debug_level)

    def make_writer(self, outfile = None, format = None):
        outfile = outfile or sys.stdout
        if format is None:
            format = self.default_output_format
        format = self.registery.normalize(format)
        return self.registery.find_writer(self, format, outfile)

    def convert(self, infile = None, outfile = None,
                input_format = None, output_format = None):
        infile = infile or sys.stdin
        outfile = outfile or sys.stdout
        
        from Martel import Dispatch
        if input_format is None:
            input_format = self.default_input_format
        input_format = self.registery.normalize(input_format)

        input_format, infile = self._get_file_format(input_format, infile)
        if input_format is None:
            raise TypeError("Could not not determine file type")
        builder = self.registery.find_builder(input_format, self)

        writer = self.make_writer(outfile, output_format)

        select_names = _get_selected_names(builder, input_format)
        parser = input_format.make_parser(select_names = select_names)

        # We can optimize StdHandler conversions
        import StdHandler
        if isinstance(builder, Dispatch.Dispatcher):
            cont_h = StdHandler.ConvertDispatchHandler(builder, writer)
        else:
            cont_h = StdHandler.ConvertHandler(builder, writer)

        writer.writeHeader()
        parser.setContentHandler(cont_h)
        parser.parseFile(infile)
        writer.writeFooter()
