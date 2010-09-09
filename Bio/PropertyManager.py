"""Stores properties associated with the class of an object (DEPRECATED).

This module is deprecated, and is expected to be removed in the next release.
If you use this module, please contact the Biopython developers via the
mailing lists.
"""
#NOTE - Adding a deprecation warning would affect Bio.Alphabet.IUPAC


# Would it be nice to have support for more than one resolver per
# class?  In the meanwhile, they could collude using a dispatch
# object.

# Do you need access to the actual resolver?

# Resolvers get the sequence because they may do a per-object lookup.

# Could cache search results for better performance.


# Dictionary which creates dictionary elements, so lookups never fail.
# The new elements are always dictionaries.
class CreateDict(dict):
    def __getitem__(self, key):
        return self.setdefault(key,{})
    
class PropertyManager:
    def __init__(self):
        self.class_property = CreateDict()
        self.class_property_resolver = CreateDict()
        self.class_resolver = {}

    def resolve(self, obj, property):
        try:
            klass = obj.__class__
        except AttributeError:
            raise KeyError("built-in instance")
        
        return self.resolve_class(klass, property)
        
    def resolve_class(self, klass, property):
        # Hopefully, we'll find the hit right away
        try:
            return self.class_property[klass][property]
        except KeyError:
            pass

        # Is there a property resolver?
        try:
            return self.class_property_resolver[klass][property](
                self, klass, property)
        except KeyError:
            pass

        # What about the class resolver?
        try:
            return self.class_resolver[klass](self, klass, property)
        except KeyError:
            pass

        # That failed, so we walk up the class tree, depth-first and
        # left-to-right (same as Python).  For each class, check if
        # the property exists, then check if the property resolver
        # exists, and finally, check for the class resolver.

        bases = list(klass.__bases__)
        while bases:
            base = bases.pop()
            try:
                return self.class_property[base][property]
            except KeyError:
                pass
            try:
                return self.class_property_resolver[base][property](
                    self, klass, property)
            except KeyError:
                pass
            try:
                return self.class_resolver[base](self, klass, property)
            except KeyError:
                pass

            # this is why the search is depth-first/right-left
            bases[:0] = list(base.__bases__)
        raise KeyError("cannot find property %s for class %s" \
                       % (property, klass))
            

default_manager = PropertyManager()
