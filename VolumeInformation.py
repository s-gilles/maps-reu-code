#!/usr/bin/python

from snappy import *

# Given the name of a manifold, returns (trace field, volume #as string) pairs in a list
# Volumes' precision is based on pari, so set it there
def get_all_volumes(m_name, engine = None):
    m = Manifold(m_name)
    try:
        v = m.ptolemy_variety(2,'all').retrieve_solutions(numerical = True)
    except: # HERE restrict to right exception type
        if engine == None:
            v = m.ptolemy_variety(2,'all').compute_solutions(numerical = True)
        else:
            v = m.ptolemy_variety(2,'all').compute_solutions(numerical = True, engine = engine)
    return list([str(g.real()) for g in [x for s in v.complex_volume_numerical() for x in s][0]]) # Structural monstrosity, sorry.

class VolumeData:

    # structure: dict poly ---> (volume, manifold)
    __init__(self, data = dict()):
        self.data = data

    def get_polys(self):
        return self.data.keys()

    def get_volumes(self,poly):
        return [rec[0] for rec in data[poly]]

    def get_manifolds(self,poly):
        return [rec[1] for rec in data[poly]]

    def get_volume_data(self,poly):
        return data[poly]

    # returns a VolumeData object containing the data from this and other; in case of a conflict, other's data takes precendence
    def combine_with(self,other):
        new_data = dict()
        new_data.update(self.data)  # HERE it occurs to me that this may do the wrong thing by overwriting records instead of combining them
        new_data.update(other.data)
        return VolumeData(data = new_data)
