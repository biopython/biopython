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

    def acquire(self, obj, prefix = ""):
        # Add the objects methods to the local tables, possibly with
        # a prefix.
        for tagname, new_method in obj._start_table.items():
            tagname = prefix + tagname
            method = _merge_methods(self._start_table.get(tagname, None),
                                    new_method, MulticallStart)
            self._start_table[tagname] = method

        for tagname, new_method in obj._end_table.items():
            tagname = prefix + tagname
            method = _merge_methods(self._end_table.get(tagname, None),
                                    new_method, MulticallEnd)
            self._end_table[tagname] = method
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


def acquire_text(dispatch, tag, attribute, verify_none = 1):
    dispatch.acquire(GetText(dispatch, attribute, verify_none), tag)
    
def acquire_append_text(dispatch, tag, attribute):
    dispatch.acquire(AppendText(dispatch, attribute), tag)
    

class Callback(DispatchHandler):
    def __init__(self, callback):
        DispatchHandler.__init__(self)
        self.callback = callback

# The 'verify_none' is to ensure that you've reset the variables
# to None before using them again.  The goal is to prevent accidental
# overwrites.  I'm not sure if this option is useful enough, and
# there is a (slight) performance hit for it.  Please me know if this
# saved you some trouble.
class GetText(DispatchHandler):
    def __init__(self, obj, attrname, verify_none = 1):
        DispatchHandler.__init__(self)
        self.obj = obj
        self.attrname = attrname
        self.verify_none = verify_none
    def start_(self, tag, attrs):  # yes, this an empty tag
        self.save_characters()     # it expects everything from the prefix
    def end_(self, tag):
        assert self.verify_none == 0 or\
               getattr(self.obj, self.attrname, None) is None, \
               "%s %r must be None (it is %s); or set verify_none to 0" % \
               (self.obj, self.attrname, getattr(self.obj, self.attrname))
        setattr(self.obj, self.attrname, self.get_characters())

class AppendText(DispatchHandler):
    def __init__(self, obj, attrname):
        DispatchHandler.__init__(self)
        self.obj = obj
        self.attrname = attrname
    def start_(self, tag, attrs):  # yes, this an empty tag
        self.save_characters()     # it expects everything from the prefix
    def end_(self, tag):
        getattr(self.obj, self.attrname).append(self.get_characters())

class SetAttr:
    def __init__(self, obj, attrname, verify_none = 1):
        self.obj = obj
        self.attrname = attrname
        self.verify_none
    def __call__(self, obj):
        assert self.verify_none == 0 or\
               getattr(self.obj, self.attrname, None) is None, \
               "%s %r must be None (it is %s); or set verify_none to 0" % \
               (self.obj, self.attrname, getattr(self.obj, self.attrname))
        setattr(self.obj, self.attrname, obj)

# Do these need weak references?
class RemapStart:
    def __init__(self, obj, new_tag):
        self.obj = obj
        self.new_tag = new_tag
    def __call__(self, tag, attrs):
        self.obj.startElement(self.new_tag, tag)
class RemapEnd:
    def __init__(self, obj, new_tag):
        self.obj = obj
        self.new_tag = new_tag
    def __call__(self, tag):
        self.obj.endElement(self.new_tag)
        

class Dispatcher(handler.ContentHandler, DispatchHandler):
    """Adapter from the standard SAX events"""
    def __init__(self, prefix = "", remap = None):
        DispatchHandler.__init__(self, prefix)

        # remap is *not* prefixed, else we couldn't remap from outside
        # our namespace.  (Is that important?)
        if remap is None:
            from expressions import Std
            remap = {"record": Std.NS}  # Is this really appropriate?
            
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
        DispatchHandler.acquire(self, obj, prefix)
        obj.setCharacterSaver(self)

    def uses_tags(self):
        # Get all the tags needed by both dictionaries
        d = {}
        d.update(self._start_table)
        d.update(self._end_table)
        return d.keys()

    # Dispatch events to the appropriate handlers
    def startElement(self, tag, attrs):
        f = self._start_table.get(tag)
        if f is not None:
            f(tag, attrs)
    def endElement(self, tag):
        f = self._end_table.get(tag)
        if f is not None:
            f(tag)

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
