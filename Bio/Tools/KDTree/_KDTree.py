# This file was created automatically by SWIG.
import _KDTreec
class KDTree:
    __setmethods__ = {}
    for _s in []: __setmethods__.update(_s.__setmethods__)
    def __setattr__(self,name,value):
        if (name == "this"):
            if isinstance(value,KDTree):
                self.__dict__[name] = value.this
                if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
                del value.thisown
                return
        method = KDTree.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value

    __getmethods__ = {}
    for _s in []: __getmethods__.update(_s.__getmethods__)
    def __getattr__(self,name):
        method = KDTree.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name

    def __init__(self,*args):
        self.this = apply(_KDTreec.new_KDTree,args)
        self.thisown = 1
    def __del__(self,_KDTreec=_KDTreec):
        if getattr(self,'thisown',0):
            _KDTreec.delete_KDTree(self)
    def set_data(*args): return apply(_KDTreec.KDTree_set_data,args)
    def build_tree(*args): return apply(_KDTreec.KDTree_build_tree,args)
    def search_center_radius(*args): return apply(_KDTreec.KDTree_search_center_radius,args)
    def get_count(*args): return apply(_KDTreec.KDTree_get_count,args)
    def neighbor_search(*args): return apply(_KDTreec.KDTree_neighbor_search,args)
    def neighbor_simple_search(*args): return apply(_KDTreec.KDTree_neighbor_simple_search,args)
    def neighbor_get_count(*args): return apply(_KDTreec.KDTree_neighbor_get_count,args)
    def get_indices(*args): return apply(_KDTreec.KDTree_get_indices,args)
    def get_radii(*args): return apply(_KDTreec.KDTree_get_radii,args)
    def neighbor_get_indices(*args): return apply(_KDTreec.KDTree_neighbor_get_indices,args)
    def neighbor_get_radii(*args): return apply(_KDTreec.KDTree_neighbor_get_radii,args)
    def __repr__(self):
        return "<C KDTree instance at %s>" % (self.this,)

class KDTreePtr(KDTree):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = KDTree
_KDTreec.KDTree_swigregister(KDTreePtr)

