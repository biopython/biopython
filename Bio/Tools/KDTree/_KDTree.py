# This file was created automatically by SWIG.
import _KDTreec
class KDTree:
    def __init__(self,*args):
        self.this = apply(_KDTreec.new_KDTree,args)
        self.thisown = 1

    def __del__(self,_KDTreec=_KDTreec):
        if self.thisown == 1 :
            _KDTreec.delete_KDTree(self)
    def set_data(*args):
        val = apply(_KDTreec.KDTree_set_data,args)
        return val
    def search_center_radius(*args):
        val = apply(_KDTreec.KDTree_search_center_radius,args)
        return val
    def get_count(*args):
        val = apply(_KDTreec.KDTree_get_count,args)
        return val
    def neighbor_search(*args):
        val = apply(_KDTreec.KDTree_neighbor_search,args)
        return val
    def neighbor_simple_search(*args):
        val = apply(_KDTreec.KDTree_neighbor_simple_search,args)
        return val
    def neighbor_get_count(*args):
        val = apply(_KDTreec.KDTree_neighbor_get_count,args)
        return val
    def get_indices(*args):
        val = apply(_KDTreec.KDTree_get_indices,args)
        return val
    def get_radii(*args):
        val = apply(_KDTreec.KDTree_get_radii,args)
        return val
    def neighbor_get_indices(*args):
        val = apply(_KDTreec.KDTree_neighbor_get_indices,args)
        return val
    def neighbor_get_radii(*args):
        val = apply(_KDTreec.KDTree_neighbor_get_radii,args)
        return val
    def __repr__(self):
        return "<C KDTree instance at %s>" % (self.this,)
class KDTreePtr(KDTree):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = KDTree





#-------------- FUNCTION WRAPPERS ------------------



#-------------- VARIABLE WRAPPERS ------------------

