#!/usr/bin/python

from snappy import *

# Given the name of a manifold, returns strings for all volumes in a list
# Set precision is based on pari, so set it there
def get_all_volumes(m_name, engine = None):
    m = Manifold(m_name)
    if engine == None:
        v = m.ptolemy_variety(2,'all').compute_solutions(numerical = True)
    else:
        v = m.ptolemy_variety(2,'all').compute_solutions(numerical = True, engine = engine)
    return list([str(g.real()) for g in [x for s in v.complex_volume_numerical() for x in s][0]]) # Structural monstrosity, sorry.
