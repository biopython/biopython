# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.
import _CKDTree
def _swig_setattr(self,class_type,name,value):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    self.__dict__[name] = value

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class KDTree(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, KDTree, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, KDTree, name)
    def __init__(self,*args):
        _swig_setattr(self, KDTree, 'this', apply(_CKDTree.new_KDTree,args))
        _swig_setattr(self, KDTree, 'thisown', 1)
    def __del__(self, destroy= _CKDTree.delete_KDTree):
        try:
            if self.thisown: destroy(self)
        except: pass
    def set_data(*args): return apply(_CKDTree.KDTree_set_data,args)
    def search_center_radius(*args): return apply(_CKDTree.KDTree_search_center_radius,args)
    def get_count(*args): return apply(_CKDTree.KDTree_get_count,args)
    def neighbor_search(*args): return apply(_CKDTree.KDTree_neighbor_search,args)
    def neighbor_simple_search(*args): return apply(_CKDTree.KDTree_neighbor_simple_search,args)
    def neighbor_get_count(*args): return apply(_CKDTree.KDTree_neighbor_get_count,args)
    def get_indices(*args): return apply(_CKDTree.KDTree_get_indices,args)
    def get_radii(*args): return apply(_CKDTree.KDTree_get_radii,args)
    def neighbor_get_indices(*args): return apply(_CKDTree.KDTree_neighbor_get_indices,args)
    def neighbor_get_radii(*args): return apply(_CKDTree.KDTree_neighbor_get_radii,args)
    def __repr__(self):
        return "<C KDTree instance at %s>" % (self.this,)

class KDTreePtr(KDTree):
    def __init__(self,this):
        _swig_setattr(self, KDTree, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, KDTree, 'thisown', 0)
        _swig_setattr(self, KDTree,self.__class__,KDTree)
_CKDTree.KDTree_swigregister(KDTreePtr)


