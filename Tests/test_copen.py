#!/usr/bin/env python

from Bio.MultiProc import copen

class T:
    def __init__(self):
        self.a = 1
        
# Can copen properly run a function with arguments?
def print_args(*args):
    for x in args:
        print x
print "opening handle"
handle = copen.copen_fn(print_args, *(range(2) + ['a', 'b', 'c']))
print "reading"
handle.read()

# Can copen properly pickle and return an object?
def obj():
    return T()
t = copen.copen_fn(obj).read()
print t.a

# Does copen properly fail when an object is not pickleable?
def unpickleable_obj():
    class NP:
        pass
    return NP()
try:
    copen.copen_fn(unpickleable_obj).read()
except AssertionError, x:
    if str(x).startswith("Error in child process"):
        print "properly failed"
    else:
        print "error"

# Test killing the process.
def take_long_time():
    import time
    time.sleep(5)
handle = copen.copen_fn(take_long_time)
print "POLL", handle.poll()   # 0
print "KILLED WITH", handle.killsig
handle.close()
print "POLL", handle.poll()   # 1
print "KILLED WITH", handle.killsig

# Test wait, poll, etc...

# Test what happens when you the function in copen_fn raises a
# traceback.
def raise_error():
    print "raising error"
    raise AssertionError, "raised"

handle = copen.copen_fn(raise_error)
try:
    handle.wait()
except AssertionError, x:
    s = str(x)   # Error in child process: ...
    print "Raised error properly"
    # chop out the path of the file, since it varies from machine to machine:
    #reported = s[:s.index('"')+1] + s[s.index('Bio/'):]
    #print reported
