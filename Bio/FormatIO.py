import sys, urllib
from xml.sax import saxutils
import ReseekFile
from Martel import Dispatch

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
    
    select_names = tuple(builder.uses_tags())
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

    def _find_builder(self, builder, resolver):
        if builder is None:
            builder = self.registery.find_builder(resolver, self)
##            if builder is None:
##                raise TypeError(
##                    "Could not find a %r builder for %r types" % \
##                    (self.name, resolver.format.name))
        return builder

    def _get_file_format(self, format, source):
        import Format        
        # By construction, this is the lowest we can go, so don't try
        if isinstance(format, Format.FormatDef):
            return format

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
        format = self.registery.normalize(format)

        format, infile = self._get_file_format(format, infile)
            
        if format is None:
            raise TypeError("Could not determine file type")

        builder = self._find_builder(builder, format)
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


    def readString(self, s, format = None, builder = None, debug_level = 0):
        if format is None:
            format = self.default_input_format
        format = self.registery.normalize(format)

        # Check to see if we need to do further identification.
        # (FormatDef is the leaf, so we don't need to go further than that.)
        import Format
        if not isinstance(format, Format.FormatDef):
            format = format.identifyString(s)
            if format is None:
                raise TypeError("Could not determine file type")

        builder = self._find_builder(builder, format)
        select_names = _get_selected_names(builder, format)
        
        if format.multirecord == 1:
            iterator = format.make_iterator("record",
                                            select_names = select_names,
                                            debug_level = debug_level)
            return FormatIOIterator(iterator.iterateString(s, builder))
        elif format.multirecord == 0:
            parser = format.make_parser(select_names = select_names,
                                        debug_level = debug_level)
            parser.setContentHandler(builder)
            parser.parseString(infile)
            return builder.document
        else:
            raise AssertionError(format.multirecord)
      
                  
    def read(self, systemID, format = None, builder = None, debug_level = 0):
        if isinstance(systemID, type("")):
            # URL
            infile = ReseekFile.ReseekFile(urllib.urlopen(systemID))
        else:
            infile = systemID
        return self.readFile(infile, format, builder, debug_level)

    def make_writer(self, outfile = sys.stdout, format = None):
        if format is None:
            format = self.default_output_format
        format = self.registery.normalize(format)
        return self.registery.find_writer(self, format, outfile)

    def convert(self, infile = sys.stdin, outfile = sys.stdout,
                input_format = None, output_format = None):
        if input_format is None:
            input_format = self.default_input_format
        input_format = self.registery.normalize(input_format)

        input_format, infile = self._get_file_format(input_format, infile)
        if input_format is None:
            raise TypeError("Could not not determine file type")

        builder = self._find_builder(None, input_format)

        writer = self.make_writer(outfile, output_format)

        select_names = _get_selected_names(builder, format)
        parser = format.make_parser(select_names = select_names)

        # We can optimize StdHandler conversions
        import StdHandler
        if isinstance(builder, Dispatch.Dispatcher):
            cont_h = StdHandler.ConvertDispatchHandler(builder, writer)
        else:
            cont_h = StdHandler.ConvertHandler(builder, writer)
        parser.setContentHandler(cont_h)
            
        parser.parseFile(infile)

