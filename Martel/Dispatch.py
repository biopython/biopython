from xml.sax import handler

# Do yourself a favor and don't use tags which mess this up.
def escape(s):
    return s.replace(":", "__")
def unescape(s):
    return s.replace("__", ":")

def _dir_class(obj):
    x = dir(obj)
    for parent in obj.__bases__:
        x = x + _dir_class(parent)
    return x


def superdir(obj):
    names = {}
    for name in dir(obj) + _dir_class(obj.__class__):
        names[name] = 1
    return names.keys()


class Multicall:
    def __init__(self, objs = None):
        if objs is None:
            self.objs = []
        else:
            self.objs = objs
    def flatten(self):
        newobjs = []
        for obj in self.objs:
            if isinstance(obj, Multicall):
                newobjs.extend(obj.flatten())
            else:
                newobjs.append(obj)
        assert newobjs, "cannot be empty"
        return newobjs
    def simplify(self):
        newobjs = self.flatten()
        if len(newobjs) == 1:
            return newobjs[0]
        return self.__class__(newobjs)
        

class MulticallStart(Multicall):
    def add(self, obj):
        # start elements are first come, first served, so append
        self.objs.append(obj)
    def __call__(self, tag, attrs):
        for obj in self.objs:
            obj(tag, attrs)

class MulticallEnd(Multicall):
    def add(self, obj):
        # end elements are last come, first served; preserves stack order
        self.objs.insert(0, obj)
    def __call__(self, tag):
        for obj in self.objs:
            obj(tag)

def _merge_methods(meth1, meth2, klass):
    if meth1 is None:
        assert meth2 is not None
        return meth2
    if meth2 is None:
        return meth1
    k = klass()
    k.add(meth1)
    k.add(meth2)
    return k.simplify()

class DispatchHandler:
    def __init__(self, prefix = ""):
        start_table = self._start_table = {}
        end_table = self._end_table = {}
        self._acquired = []
        self._prefix = prefix
        self.supported_features = []

        for methodname in superdir(self):
            if methodname.startswith("start_"):
                method = getattr(self, methodname)
                escaped_tagname = methodname[6:]
                tagname = prefix + unescape(escaped_tagname)
                assert not start_table.has_key(tagname)  # by construction
                start_table[tagname] = method
            elif methodname.startswith("end_"):
                method = getattr(self, methodname)
                escaped_tagname = methodname[4:]
                tagname = prefix + unescape(escaped_tagname)
                assert not end_table.has_key(tagname)  # by construction
                end_table[tagname] = method

    def get_supported_features(self):
        return self.supported_features
        

    def acquire(self, obj, prefix = ""):
        # Add the objects methods to the local tables, possibly with
        # a prefix.
        for tagname, new_method in obj._start_table.items():
            tagname = prefix + tagname

            # If you get an AttributeError exception with the text
            #     "Show instance has no attribute '_start_table'"
            # then you forgot to initialize the parent class
            method = _merge_methods(self._start_table.get(tagname, None),
                                    new_method, MulticallStart)
            
            self._start_table[tagname] = method

        for tagname, new_method in obj._end_table.items():
            tagname = prefix + tagname
            method = _merge_methods(self._end_table.get(tagname, None),
                                    new_method, MulticallEnd)
            self._end_table[tagname] = method

        # This isn't technically correct because one acquisition might
        # support a feature while another might not, so in that case
        # we need disable feature support.  Getting that information
        # is tricky, so I decided that for now if you want feature
        # support you better know what you're doing.  :(
        d = {}
        for x in self.supported_features + obj.get_supported_features():
            d[x] = 1
        self.supported_features[:] = d.keys()
        self._acquired.append(obj)

    def setCharacterSaver(self, saver):
        # Is there a cycle here?  use a weak object?
        self.save_characters = saver.save_characters
        saver.save_characters()
        self.get_characters = saver.get_characters
        saver.get_characters()
        for obj in self._acquired:
            obj.setCharacterSaver(saver)
            
    def save_characters(self):
        raise AssertionError("Not yet set by a Dispatcher")

    def get_characters(self):
        raise AssertionError("Not yet set by a Dispatcher")


class Callback(DispatchHandler):
    def __init__(self, callback):
        DispatchHandler.__init__(self)
        self.callback = callback

# Do these need weak references?
class RemapStart:
    def __init__(self, obj, new_tag):
        self.obj = obj
        self.new_tag = new_tag
    def __call__(self, tag, attrs):
        self.obj.startElement(self.new_tag, attrs)
class RemapEnd:
    def __init__(self, obj, new_tag):
        self.obj = obj
        self.new_tag = new_tag
    def __call__(self, tag):
        self.obj.endElement(self.new_tag)
        

class Dispatcher(handler.ContentHandler, DispatchHandler):
    """Adapter from the standard SAX events"""
    def __init__(self, prefix = "", remap = {}):
        DispatchHandler.__init__(self, prefix)

        # remap is *not* prefixed, else we couldn't remap from outside
        # our namespace.  (Is that important?)
        for old_tagname, new_tagname in remap.items():
            start = RemapStart(self, new_tagname)
            method = _merge_methods(self._start_table.get(old_tagname, None),
                                    start, MulticallStart)
            self._start_table[old_tagname] = method

            end = RemapEnd(self, new_tagname)
            method = _merge_methods(self._end_table.get(old_tagname, None),
                                    end, MulticallEnd)
            self._end_table[old_tagname] = method
                                    
            
        self._save_stack = []
        self.setCharacterSaver(self)

    def acquire(self, obj, prefix = ""):
        if prefix[:1] == ":":
            prefix = prefix[1:]
        else:
            prefix = self._prefix + prefix
        DispatchHandler.acquire(self, obj, prefix)
        obj.setCharacterSaver(self)

    def uses_tags(self):
        # Get all the tags needed by both dictionaries
        d = {}
        d.update(self._start_table)
        d.update(self._end_table)
        return d.keys()

    # Dispatch events to the appropriate handlers
    def startDocument(self):
        self.startElement(self._prefix, {})
            
    def startElement(self, tag, attrs):
        f = self._start_table.get(tag)
        if f is not None:
            f(tag, attrs)
    def endElement(self, tag):
        f = self._end_table.get(tag)
        if f is not None:
            f(tag)

    def endDocument(self):
        self.endElement(self._prefix)

    # Stack-based way for the handlers to get data without
    # stepping on each other's toes -- so long as the calls
    # are balanced!
    def save_characters(self):
        if self._save_stack:
            self._save_stack.append(len(self._save_text))
        else:
            self._save_text = ""
            self._save_stack.append(0)
    def get_characters(self):
        pos = self._save_stack.pop()
        return self._save_text[pos:]

    def characters(self, s):
        if self._save_stack:
            self._save_text += s
