from Bio.PropertyManager import PropertyManager

import sys
if sys.version_info[0] >= 3:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(\
        "This deprecated module doesn't work on Python 3.")

def test():
    pm = PropertyManager()
    class Foo:
        pass
    class Bar:
        pass
    class FooBar(Foo, Bar):
        pass
 
    data = [4, 5, 6]
    pm.class_property[Foo]["name"] = "only Foo"
    pm.class_property[Bar]["time"] = "11 pm"
    pm.class_property[FooBar]["list"] = data
    pm.class_property[FooBar]["name"] = "Foo Bar"
 
    f = Foo()
    b = Bar()
    fb = FooBar()
    print pm.resolve(f, "name")
    #print pm.resolve(b, "name")
    print pm.resolve(fb, "name")
 
    def list_resolver(manager, klass, property, FooBar = FooBar):
        print "resolving list"
        x = manager.resolve_class(FooBar, "list")
        y = []
        for a in x:
            y.append(a + 10)
        return y
 
    pm.class_property_resolver[Bar]["list"] = list_resolver
    print pm.resolve(fb, "list")
    print pm.resolve(b, "list")
    print pm.resolve(b, "list")
    data.append(-10)
    print pm.resolve(b, "list")
 
    def prop_resolver(manager, klass, property):
        print "I am called"
        x = str(klass) + property
        manager.class_property[klass][property] = x
        return x
 
    pm.class_resolver[Foo] = prop_resolver
    print pm.resolve(f, "qwq")
    print pm.resolve(f, "qwq")

test()
