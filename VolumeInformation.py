from snappy import *

def get_all_volumes(m_name, engine = None):
    """
    Given the name of a manifold (as a string), returns the a list of
    [the real part of] each complex volume that is a solution of the
    ptolemy equations for the manifold. The volumes are all returned as
    strings.
    """
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
    """
    A structure that contains volumes and their accompanying manifolds
    for each invariant trace field polynomial
    """
    # structure: dict poly ---> (volume, manifold)
    def __init__(self, data = dict()):
        self.data = data

    def get_polys(self):
        """
        Returns a list of all polynomials with data held by this
        VolumeData
        """
        return self.data.keys()

    def get_volumes(self,poly):
        """
        Returns a list of all volumes for this polynomial known to this
        VolumeData
        """
        return [rec[0] for rec in self.data[poly]]

    def get_manifolds(self,poly):
        """
        Returns a list of all manifolds for this polynomial known to
        this VolumeData
        """
        return [rec[1] for rec in self.data[poly]]

    def get_volume_data(self,poly):
        """
        Returns a list containing tuples, each of which is (v,m),
        representing a volume and a manifold. This list contains all
        such (volume, manifold) pairings for this polynomial known to
        this VolumeData
        """
        return self.data[poly]

    # returns a VolumeData object containing the data from this and other; in case of a conflict, other's data takes precendence
    def combine_with(self,other):
        """
        Returns a new VolumeData object. This VolumeData object contains
        all (volume, manifold) pairs (and their associated polynomials)
        which are contained in either this VolumeData object or other.

        Note: if this VolumeData and other contain contradictory
        information, both will be stored in the resultant VolumeData
        object. This may result in, for example, two different tuples
        recording slightly different volumes for the same manifold under
        the same polynomial.
        """
        new_data = dict()

        all_polys = set(self.data.keys() + other.data.keys())
        for poly in all_polys:
            new_data[poly] = self.data.get(poly, list()) + other.data.get(poly, list())
        return VolumeData(data = new_data)
