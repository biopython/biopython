import weakref, os, sys
import Format, _FmtUtils

def _issubtype(child, parent):
    if child.name == parent.name:
        return 1
    for term in child._parents:
        if _issubtype(term, parent):
            return 1
    return 0


def _find_child(children, name, parent_name, child_name, direction):
    i = 0
    for child in children:
        if child[0].name == name:
            return i
        i = i + 1
    raise TypeError("%r not defined in %r, so cannot put %r %s it" % \
                    (name, parent_name, child_name, direction))


class FormatRegistry:
    def __init__(self, loadpath):
        self._name_table = {}
        self._abbrev_table = {}
        self.loadpath = loadpath
        self._autoloaded = 0
        self._autoloading = 0

    def _autoload(self):
        if self._autoloaded or self._autoloading:
            return
        self._autoloading = 1
        _FmtUtils.load_basemodule(self.loadpath)
        self.autoloading = 0
        self._autoloaded = 1
        
    def register_format(self, **kwargs):
        self._autoload()
        if kwargs.has_key("expression"):
            format = Format.FormatDef(**kwargs)
        else:
            format = Format.FormatGroup(**kwargs)
        name = format.name
        abbrev = format.abbrev
        if self._name_table.has_key(name):
            raise TypeError("%r is a duplicate entry" % (name,))
        if self._abbrev_table.has_key(abbrev):
            raise TypeError("%r is a duplicate entry" % (abbrev,))
        
        self._name_table[name] = format
        self._abbrev_table[abbrev] = format

    def __getitem__(self, name):
        self._autoload()
        return self._name_table[name]  # raises KeyError for unknown formats

    def get(self, name, default = None):
        self._autoload()
        return self._name_table.get(name, default)

    def keys(self):
        self._autoload()
        return self._name_table.keys()
    def values(self):
        self._autoload()
        return self._name_table.values()
    def items(self):
        self._autoload()
        return self._name_table.items()
    

    def normalize(self, name_or_format): # XXX appropriate?
        if isinstance(name_or_format, type("")):
            # It's a name
            return self[name_or_format]
        return name_or_format

    def link(self, parent, child, filter = None, before = None, after = None,
             parent_before = None, parent_after = None):
        self._autoload()
        if not self._name_table.has_key(parent):
            raise TypeError("parent (%s) not found in registery" % (parent,))
        parent = self._name_table[parent]
        
        if not self._name_table.has_key(child):
            raise TypeError("child (%s) not found in registery" % (child,))
        child = self._name_table[child]

        if _issubtype(parent, child):
            raise TypeError("will not add a cycle between %r and %r" % \
                            (parent.name, child.name))
        
        if isinstance(filter, type("")):
            filter = _FmtUtils.CacheParser(_FmtUtils.ExpressionLoader(filter))

        children = parent._children
        new_term = (weakref.proxy(child), filter)

        if before is None:
            if after is None:
                # neither defined, so stick at the end
                children.append(new_term)
            elif after == "*end*":
                children.append(new_term)
            else:
                pos = _find_child(children, after,
                                  parent.name, child.name,
                                  "after") + 1
                children.insert(pos, new_term)
        else:
            if after is not None:
                raise TypeError("Cannot specify both 'before' and 'after'")
            elif before == "*start*":
                children.insert(0, new_term)
            else:
                pos = _find_child(children, before,
                                  parent.name, child.name,
                                  "before")
                children.insert(pos, new_term)
        

        parents = child._parents
        if parent_before is None:
            if parent_after is None:
                parents.append(parent)
            elif parent_after == "*end*":
                parents.append(parent)
            else:
                pos = _find_child(parents, parent_after,
                                  child.name, parent.name, "parent_after") + 1
                parents.insert(pos, child)
        else:
            if parent_after is not None:
                raise TypeError(
                    "Cannot specify both 'parent_before' and 'parent_after'")
            elif before == "*start*":
                parents.insert(0, child)
            else:
                pos = _find_child(parents, parent_before,
                                  child.name, parent.name, "parent_before")
                parents.insert(pos, child)
    
    def find_builder(self, from_format, to_io):
        from_format_name = from_format.abbrev
        to_class_name = to_io.abbrev
        
        modulename = self.loadpath + ".builders." + to_class_name
        module = _FmtUtils.load_module(modulename)

        name = _FmtUtils._load_possible_module(modulename,
                            from_format._get_parents_in_depth_order())
        if name is None:
            raise TypeError("Cannot find builder for %r" % \
                            (to_class_name,))
        module = getattr(module, name)
        
        return module.make_builder()  # Get the builder function

    def find_writer(self, from_io, to_format, outfile):
        from_class_name = from_io.abbrev

        modulename = self.loadpath + ".writers." + from_class_name
        module = _FmtUtils.load_module(modulename)

        name = _FmtUtils._load_possible_module(modulename,
                               to_format._get_children_in_depth_order())
        if name is None:
            raise TypeError("Cannot find writer for %r" % \
                            (from_class_name,))
        module = getattr(module, name)

        return module.make_writer(outfile)
