# Various helper routines for the Bioformat IO system.

# I don't think any of these routines are generally useful but I
# didn't know where else to put them.

import string, sys, os

def load_module(modulename):
    try:
        module = __import__(modulename)
    except SyntaxError, exc:
##        raise SyntaxError("%s during import of %r" % (exc, modulename)), \
##              None, sys.exc_info()[2]
        raise
    except ImportError, exc:
        raise ImportError("%s during import of %r" % (exc, modulename)), \
              None, sys.exc_info()[2]
    for name in modulename.split(".")[1:]:
        module = getattr(module, name)
    return module

def load_object(path):
    terms = string.split(path, ".")
    s = terms[0]
    # Import all the needed modules
    # (Don't know which are modules and which are classes, so simply
    # stop when imports fail.)
    # The order of appends is correct, since the last element cannot
    # be a module.
    x = __import__(s)
    prev_term = s
    for term in terms[1:]:
        try:
            __import__(s)
        except SyntaxError, exc:
##            raise SyntaxError("%s during import of %r" % (exc, s)), \
##                  None, sys.exc_info()[2]
            raise
        except ImportError, exc:
            # This is the only way I know to tell if the module
            # could not be loaded because it doesn't exist.
            error_text = str(exc)
            if error_text.find("No module named %s" % prev_term) == -1:
                raise
            break
        if not term:
            raise TypeError("There's a '.' in the wrong place: %r" % \
                            (path,))
        s = s + "." + term
        prev_term = term

    # Get the requested object
    s = terms[0]
    for term in terms[1:]:
        try:
            x = getattr(x, term)
        except AttributeError:
            raise AttributeError("%s object (%r) has no attribute %r" % \
                                 (type(x).__name__, s, term))
        s = s + "." + term
    return x

def load_expression(path):
    from Martel import Expression
    x = load_object(path)
    if x is not None:
        if not isinstance(x, Expression.Expression):
            try:
                klass = x.__class__.__name__
            except AttributeError:
                klass = type(x)
                raise TypeError("%r should be a Martel Expression but " \
                                "is a %r" % (path, klass))
        return x
 
    # Expression not found; make a useful error message
    msg = "Could not find %r\n" % (path,)
    msg = msg + "(You may need to add the top-level module to the PYTHONPATH)"
    raise TypeError(msg)

def load_basemodule(modulename):
    modulename = modulename + ".formatdefs"
    module = load_module(modulename)
    filename = module.__file__
    dirname = os.path.dirname(filename)
    filenames = []
    for filename in os.listdir(dirname):
        if filename[:1] != "_" and \
           (filename[-3:] == ".py" or filename[-4:] in (".pyc", ".pyo")):
            filenames.append(filename)
    for filename in filenames:
        import_name = modulename + "." + os.path.splitext(filename)[0]
        # Can do multiple loads, but that's okay since it's cached
        # in sys.modules
        load_module(import_name)

def _load_possible_module(basemodulename, possible_formats):
    for format in possible_formats:
        try:
            __import__(basemodulename + "." + format.abbrev)
        except ImportError, exc:
            # This is the only way I know to tell if the module
            # could not be loaded because it doesn't exist.
            error_text = str(exc)
            if error_text.find("No module named %s" % format.abbrev) == -1:
                raise
            continue
        return format.name
    return None


class ExpressionLoader:
    def __init__(self, expression_path):
        self.expression_path = expression_path

    def __getattr__(self, name):
        if name[:2] == "__":
            raise AttributeError(name)
        expression = load_expression(self.expression_path)

        # Perform a bit of Python magic and become the other object
        self.__dict__.update(expression.__dict__)
        self.__class__ = expression.__class__
        return getattr(self, name)
        

# Only caches parsers for make_parser, not iterators
class CacheParser:
    def __init__(self, expression):
        self.expression = expression
        self._parsers = {}

    def make_parser(self, debug_level = 0):
        if self._parsers.get(debug_level) is None:
            parser = self.expression.make_parser(debug_level = debug_level)
            self._parsers[debug_level] = parser
        return self._parsers[debug_level].copy()
