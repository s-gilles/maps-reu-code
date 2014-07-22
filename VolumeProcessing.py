#!/usr/bin/python

import copy
import re
import sys
import traceback
import fractions

from numpy.linalg import det
from SpanFinder import find_span
from PseudoVols import VolumeData, get_potential_trace_fields
from cypari import *
from fractions import Fraction
from itertools import combinations
from snappy import *

EPSILON = .0000000000001
MAX_COEFF = 4096
MAX_ITF = 8
SOL_TYPE_STRINGS = ['not_attempted', 'geometric', 'nongeometric', 'flat', 'degenerate', 'unrecognized', 'none_found']   # globalize
gen.pari.set_real_precision(100)

# This class is just a wrapper for the structure storing polynomial/volume data.
# Having it avoids opaque references to the particular way data is stored that might change in the future.

class Dataset:
    """
    A class representing a collection of computed volumes. This is
    really a wrapper around a nested structure of dictionaries of
    dictionaries of.

    This class should probably not be constructed directly.  Rather, it
    should be created through a method like read_csv()
    """
    def __init__(self, data_dict = dict()):
        self.data = data_dict

    def __str__(self):
        return str(self.data)

    def copy(self):
        """
        Return another dataset containing a copy of the data that this
        Dataset represents.
        """
        return Dataset(copy.deepcopy(self.data))

    def get_polys(self):
        """
        Return all trace fields as polynomials (in string form) recorded
        by this Dataset.
        """
        return self.data.keys()

    def get_roots(self,poly):
        """
        Return all roots of a given trace field (as a polynomial in
        string form) recorded by this Dataset.

        Raises KeyError if the polynomial is not in the Dataset.
        """
        return self.data[poly][0].keys()

    def get_degree(self,poly):
        """Get the degree of this polynomial (input as a string) inside this
        Dataset.  This should be equivalent to

        gen.pari(poly).poldegree()

        Discrepancies may reveal issues with the data that constructed
        this Dataset.
        """
        return self.data[poly][1]

    def get_ncp(self,poly):
        """
        Get the number of complex places of this polynomial (input as a
        string) recorded by this Dataset.
        """
        return self.data[poly][2]

    def get_disc(self,poly):
        """
        Get the discriminant of this polynomial (input as a string)
        recorded by this Dataset. This should be equal to

        gen.pari(poly).nfdisc()

        Discrepancies may reveal issues with the data that
        constructed this Dataset.
        """
        return self.data[poly][3]

    def get_factored_disc(self,poly):
        """
        Get the factored discriminant of this polynomial, as a
        human-readable string of prime power factorization. The string
        it returns should (in the obvious way) represent the same
        factorization as

        gen.pari(pol).nfdisc().factor()

        Discrepancies may reveal issues with the data that
        constructed this Dataset.
        """
        return self.data[poly][4]

    def get_volumes(self,poly,root):
        """
        Return, as a list of strings, each volume for the given
        polynomial (input as string) and root (input as string) recorded
        in this Dataset.
        """
        return self.data[poly][0][root].keys()

    def has_poly(self,poly):
        """
        Returns true if there is any volume stored for the given
        polynomial (input as a string) in this Dataset.
        """
        return poly in self.data.keys()

    # returns a geometric manifold's record, or None failing that
    def get_geom_manifold(self,poly,root,vol):
        """
        Given a polynomial, a root, and a volume (as strings), return either a triple
        (manifold name, number of simplices, solution type) where

        manifold name is a string

        number of simplices is an integer reflecting the triangulation
        that was recorded in this Dataset

        solution type is the string 'geometric')

        If no such triple can be found, (None, None, None) is returned.
        """
        for rec in self.get_manifold_data(poly,root,vol):
            if rec[2] == 'geometric':
                return rec
        return (None,None,None)

    def get_manifold_data(self,poly,root,vol):
        """
        Returns an arbitrarily chosen triple (manifold name, number of
        simplices, solution type) for a given polynomial, root, and
        volume, where

        manifold name is a string

        number of simplices is an integer reflecting the triangulation
        that was recorded in this Dataset

        solution type is a string corresponding to a short description
        of Manifold.get_solution_type() for the triangulation that was
        recorded in this Dataset.
        """
        return self.data[poly][0][root][vol][0]

    def get_pared_manifolds(self,poly,root,vol):
        """
        When a Dataset is pared, multiple manifolds that meet the same
        polynomial, root, volume triple are compressed, with only one
        representative volume stored.  The names of manifolds which were
        not chosen to be this representative volume are also recorded,
        and may be retrieved by this method.
        """
        return self.data[poly][0][root][vol][1]

    def remove_volume(self,poly,root,vol):
        """
        Remove a volume from the dataset.  This will not work if the
        volume has been already removed by paring. If the last volume
        for a root, or the last root for a polynomial is removed, the
        higher-level element will be removed as well.
        """
        rec = self.data[poly][0][root].get(vol)
        del self.data[poly][0][root][vol]
        return rec

    def pare_all_volumes(self):
        """
        Compress all matching (polynomial, root, volume) triples, so
        that only one representative manifold is stored.  Manifolds
        which are discarded have their names stored, and can be
        retrieved by get_pared_volumes()
        """
        for p in self.get_polys():
            for r in self.get_roots(p):
                for v in self.get_volumes(p,r):
                    self.pare_volume(data,p,r,v)

    def pare_volume(self,poly,root,vol):
        """
        Compress a (polynomial, root, volume) triple, so that only one
        representative manifold is stored for it.  Manifolds which are
        discarded have their names stored, and can be retrieved by
        get_pared_volumes()
        """
        mdata = self..get_manifold_data(poly,root,vol)
        mpared = self.get_pared_manifolds(poly,root,vol)
        while len(mdata) > 1:
            mpared.append(mdata.pop(1)[0])

    def cull_all_volumes(self, epsilon = EPSILON):
        """
        Remove all volumes that are integer multiples of another,
        smaller volume. These volumes are not pared (and so cannot be
        retrieved by get_pared_volumes()), they are removed outright
        from the Dataset. For large Datasests, this may free resources
        and make dealing with the Dataset faster.
        """
        for p in self.get_polys():
            for r in self.get_roots(p):
                self.cull_volumes(p,r,epsilon = epsilon)

    def cull_volumes(self,poly,root,epsilon = EPSILON):
        vols = self.get_volumes(poly,root)
        # vols = self.self[poly][0][root].keys()
        i = 0
        while i < len(vols) - 1:
            j = i + 1
            while j < len(vols):
                try:
                    if is_int(float(vols[i])/float(vols[j]), epsilon = epsilon) and gen.pari(vols[i] + ' > ' + vols[j]):
                        # We have to throw away (culled) manifold names to let all culled manifolds have the same volume
                        # [j] divides [i] so remove [i]
                        self.remove_volume(poly,root,vols.pop(i))
                        # i is already effectivley incremented, so we must offset it
                        i = i-1
                        break
                    elif is_int(float(vols[j])/float(vols[i]), epsilon = epsilon):
                        # this time, remove [j]
                        self.remove_volume(poly,root,vols.pop(j))
                        # j is effectivley incremented, no need to do it
                    else:
                        j += 1
                    except (ValueError, ZeroDivisionError): # bad quotient; not a linear combination either way so...
                    j += 1
            i += 1



    # Returns a Dataset with the data from self and other; in case of a conflict, other's values beat self's or both are kept
    # Therefore, one is advised to use this on disjoint Datasets or pare volumes afterwards
    def combine_with(self,other):
        new_data = dict(self.data)
        for p in other.get_polys():
            new_data.setdefault(p,[dict(),other.get_degree(p),other.get_ncp(p),other.get_disc(p),other.get_factored_disc(p)])
            for r in other.get_roots(p):
                new_data[p][0].setdefault(r,dict())
                for v in other.get_volumes(p,r):
                    new_data[p][0][r].setdefault(v,[list(),list()])
                    new_data[p][0][r][v][0].extend(other.get_manifold_data(p,r,v))
                    new_data[p][0][r][v][1].extend(other.get_pared_manifolds(p,r,v))
                    for dim in new_data[p][0][r][v]:
                        dim = list(set(dim))    # Remove duplicates.
        return Dataset(data_dict = new_data)

    def get_representative_element(self, poly, root):
        """
        Return a triplet containing data for the manifold of smallest
        volume with the given field
        """
        minvol = (None, sys.float_info.max)
        for v in self.get_volumes(poly,root):
            try:
                if float(v) < minvol[1]:
                    minvol = (v, float(v))
            except ValueError:
                continue    # probably v=0+; doesn't really matter which we pick anyway
        if minvol[0] is None:
            return None # means no geometric solutions were found for this field
        else:
            for m in self.get_manifold_data(poly,root,minvol[0]):
                return (minvol[0],m)

    def get_representative_dataset(self):
        """
        Returns a new Dataset,
        """
        newdata = dict()
        for poly in self.get_polys():
            newdata[poly] = [dict()]+self.data[poly][1:]    # initialize list of volumes to be empty
            for root in self.get_roots(poly):
                md = self.get_representative_element(poly,root)
                if md is not None:  # we actually have something geometric for this root
                    newdata[poly][0][root] = {md[0] : [[md[1]],list()]}
            if newdata[poly][0] == dict():  # no roots gave us geometric solutions
                del newdata[poly]
        return Dataset(data_dict = newdata)

    # Returns false if contents look very wrong (no x in polynomial slot, etc.)
    # Only checks very shallowly for the first record data.__iter__.next() returns, so no guarantees
    def sane(self):
        """
        Returns false if contents look very wrong (no x in polynomial
        slot, etc.). This is only a shallow check, and this function
        returning True is no guarantee that the Dataset has no errors.
        """
        try:
            try:
                p = self.data.keys().__iter__().next()
            except StopIteration:
                return True # empty Dataset is sane
            if 'x' not in p:
                return False
            r = self.get_roots(p).__iter__().next()
            if 'I' not in r and 'i' not in r:
                return False
            try:
                if not is_int(float(self.get_degree(p))):
                    return False
                if not is_int(float(self.get_ncp(p))):
                    return False
                if not is_int(float(self.get_disc(p))):
                    return False
                if '1' not in self.get_factored_disc(p) and '^' not in self.get_factored_disc(p):
                    return False
                v = self.get_volumes(p,r).__iter__().next()
                float(v)    # for the below
            except ValueError:
                return False    # something should have been a float and wasn't
            m = self.get_manifold_data(p,r,v)
            if m[0][0][-1] == ',':
                return False    # probably a Dehn surgery got spliced
            if not is_int(float(m[0][1])):
                return False
            # testing solution type would be annoying
        except:
            return False    # unexpected errors probably mean we aren't sane
        return True

    def remove_non_geometric_elements(self):
        """
        Delete all volumes if they have any solution type except
        geometric.
        """
        to_kill = list()
        for poly,polyinf in self.data.items():
            polydict = polyinf[0]
            for root,rootdict in polydict.items():
                for vol,manifolds in rootdict.items():
                    manifolds[0] = [ m for m in manifolds[0] if m[2] == 'geometric' ]
                    if not manifolds[0]:
                        to_kill.append((poly,root,vol))
        for poly,root,vol in to_kill:
            self.remove_volume(poly,root,vol)
        for p in self.get_polys():  # These blank parts probably do no harm, but just to be sure...
            for r in self.get_roots(p):
                if not self.get_volumes(p,r):
                    del self.data[p][0][r]
            if not self.get_roots(p):
                del self.data[p]

    def get_nice_manifold_name(self,poly,root,vol):
        """
        For a (polynomial, root, volume) triple (all input a strings),
        return the ``nicest'' matching manifold.
        """
        nms = [rec[0] for rec in self.get_manifold_data(poly,root,vol)]
        nms.extend(self.get_pared_manifolds(poly,root,vol))
        opt = ('', sys.float_info.max)
        for nm in nms:
            if 0 < len(nm) < opt[1]:    # TODO: finish _niceness and replace len with it; why is 0 < required?
                opt = (nm,len(nm))
        return opt[0]

    def smush_volumes(self, epsilon = EPSILON)
        """
        If some volumes differ by less than epsilon, combine them,
        keeping one of them arbitrarily.

        Note: This method is slightly nondeterministic; may not
        get large, dense (compared to epsilon) clumps if iteration order is
        unfavorable.  However, that is very unlikely to happen, and when it is
        (lattice generators very close to each other) you don't want smush
        to work anyway.

        Note: This method is supreceded by cull_all_volumes(), and has
        only been briefly tested.
        """:
        d = self.copy()
        balls = list()
        for p in d.get_polys():
            for r in d.get_roots(p):
                vol_data = d.data[p][0][r]
                balls = list()
                vols = list(vol_data.keys())
                for v in vols:  # find close volumes
                    for w in [w for w in vols if w is not v]:
                        if abs(float(w)-float(v)) < epsilon:
                            balls.append(set([v,w]))
                def _br(balls): # combine balls by finding graph components
                    for b in balls:
                        for u in [u for u in balls if u is not b]:
                            if not b.isdisjoint(u):
                                balls.remove(b)
                                balls.remove(u)
                                balls.append(b.union(u))
                                return True # acted, might have more to do; to avoid mutation errors...
                    return False
                while _br(balls):           # ...do this
                    pass
                for b in balls: # combine data for each volume into one
                    nrec = [list(),list()]
                    n = 0
                    for v in b:
                        nrec[0].extend(vol_data[v][0])  # Hopefully there should be no duplicates in this case.
                        nrec[1].extend(vol_data[v][1])
                        del vol_data[v]
                    vol_data[v] = nrec  # bit of an abuse of v
        return d

    # Given a valid Manifold object or manifold name, returns whatever we have on it
    # Form: [InvariantTraceField,Root,Volume,SolutionType,GeomAlternative,NiceAlternative] or None if it couldn't be found
    def search_for_manifold(self,man):
        man = str(man)
        for p in self.get_polys():
            for r in self.get_roots(p):
                for v in self.get_volumes(p,r):
                    if man in [rec[0] for rec in self.get_manifold_data(p,r,v)] or man in self.get_pared_manifolds(p,r,v):
                        out = [p,r,v,None,None,None]
                        if man in [rec[0] for rec in self.get_manifold_data(p,r,v)]: # search ourselves (save vs. randomize())
                            for rec in self.get_manifold_data(p,r,v):
                                if man in rec:
                                    out[3] = rec[2]
                        else:
                            out[3] = SOL_TYPE_STRINGS[int(Manifold(man).solution_type(enum = True))]
                        if out[3] == 'geometric':
                            out[4] = man
                        else:
                            out[4] = self.get_geom_manifold(p,r,v)
                        out[5] = self.get_nice_manifold_name(p,r,v)
                        return out
        return None

def _niceness(nm):
    n = 0.0
    k = nm[0]
    if k == 'm':
        n += 0
    elif k == 'v':
        n += 1
    elif k == 't':
        n += 2
    elif k == 'o':
        n += 3
    elif k == 'K':  # TODO: apply knot and link # penalties
        n += 4
    elif k == 'L':
        n += 5
    elif k == 'b':
        n += 6
    elif k == 'B':
        n += 7
    elif k == 'D':
        n += 8
    else:
        n += 10
    # n *= gap TODO
    # apply dehn penalties TODO
    return n

def quick_read_csv(filenm, seperator = ';', sub_seperator = '|'):
    try:
        f = open(filenm,'r')
        f.readline()
        d = read_csv(f, seperator = seperator, sub_seperator = sub_seperator)
        f.close()
        return d
    except:
        if f:
            f.close()
        raise


# combines two output files from this program
def quick_combine_files(filenms, fileseps, out_filenm, out_seperator = ';', out_append = False):
    dsets = list()
    for i in xrange(len(filenms)):
        inf = open(filenms[i],'r')
        try:
            inf.readline()  # skip header
            dsets.append(read_csv(inf, seperator = fileseps[i]))
        finally:
            inf.close()
    for d in dsets[1:]:
        dsets[0].combine_with(d)
    if out_append:
        ouf = open(out_filenm, 'a')
    else:
        ouf = open(out_filenm, 'w')
    write_csv(ouf, dsets[0], seperator = out_seperator, append = out_append)

# read in raw csv in_file, pare and cull it, write it to out_file
def quick_preprocess(in_filenm, out_filenm, in_seperator = ';', out_seperator = ';', out_append = False):
    inf = open(in_filenm,'r')
    try:
        inf.readline()  # skip header
        d = read_raw_csv(inf, seperator = in_seperator)
    finally:
        inf.close()
    pare_all_volumes(d)
    cull_all_volumes(d)
    if out_append:
        ouf = open(out_filenm,'a')
    else:
        ouf = open(out_filenm,'w')
    try:
        write_csv(ouf, d, seperator = out_seperator, append = out_append)
    finally:
        ouf.close()

# Load a CSV file organized by manifold and reorganize it by polynomial and volume.
# The result: dict poly ---> (dict roots ----> (dict vols ---> (list (manifold name, tetrahedra, soltype), list (pared manifolds)), degree etc.)
def read_raw_csv_from_file(in_file, seperator = ';'):
    return read_raw_csv(in_file.readlines(), seperator)

# Returns true if the given strings are equal or complex conjugates as formatted by snap: a+b*I, a-b*I
def _up_to_conjugates(z,w):
    zp = re.findall(r'([+-]?[\d.]+)',z)
    wp = re.findall(r'([+-]?[\d.]+)',w)
    return len(zp) == len(wp) == 2 and zp[0] == wp[0] and up_to_sign(zp[1],wp[1])

# Returns true if one of the strings is just -the other
# Should only be applied to non sci-notation floats' strings
def up_to_sign(x,y):
    return re.search(r'[\d.]+',x).group() == re.search(r'[\d.]+',y).group()

# Given a+b*I, returns a\pm b*I
def _get_conjs(z):
    return z[:1]+z[1:].replace('+','\xb1').replace('-','\xb1')

def read_raw_csv(contents, seperator = ';'):
    data = dict()
    # Obviously this code is highly sensative to any changes in the output format of VolumeFinder.py
    for l in contents:
        l = l.replace(' ', '')

        if seperator == ',':    # special cased since ',' appears in Dehn surgery
            w = re.findall('"([^"]*)"', l)
        else:
            w = l.replace('\n','').replace('"','').split(seperator)

        # Since order got changed (for some unknown reason):
        try:
            w = [w[0],w[9],w[4],w[1],w[5],w[2],w[6],w[3],w[7],w[8]]
        except:
            print('Error with line ' + str(l))
            continue
        # Incase the disc was 1, a temporary hack:
        # if len(w) == 8:
        #   w.append('')
        # w[0]: manifold name ---------------------------> m[0] for m in data[poly][0][root][vol][0]
        # w[1]: manifold simplices ----------------------> m[1] for m in data[poly][0][root][vol][0]
        # w[2]: volume ----------------------------------> v in data[poly][0][root].keys()
        # w[3]: invariant trace field polynomial --------> p in data.keys()
        # w[4]: polynomial degree -----------------------> data[poly][1]
        # w[5]: polynomial root -------------------------> r in data[poly][0].keys()
        # w[6]: manifold solution type ------------------> m[2] for m in data[poly][0][root][vol][0]
        # w[7]: polynomial number of complex places -----> data[poly][2]
        # w[8]: polynomial discriminant -----------------> data[poly][3]
        # w[9]: polynomial discriminant (factorized) ----> data[poly][4]
        # vr = data.setdefault(w[3],[dict(),w[4]])[0].setdefault(w[2],[list(),list(),w[5]])[0].append(w[0:2]) # OLD
        # # why was vr set just now and not used?
        vol_entry = data.setdefault(w[3],[dict(),w[4]])[0].setdefault(w[5],dict()).setdefault(w[2],[list(),list()])

        vol_entry[0].append((w[0],w[1],w[6]))

        if len(data[w[3]]) == 2:
            data[w[3]].extend(w[7:10])
    return Dataset(data)

# Reads a CSV produced by write_csv and returns the contents as a dataset object
# This variant handles csv's before we swapped column order around a bit.
def read_old_csv(in_file, seperator = ';'):
    data = dict()
    for l in in_file.readlines():
        if seperator == ',':    # again special cased
            w = re.findall('"([^"]*)"', l)
        else:
            w = l.replace('\n','').replace('"','').split(seperator)
        vol_entry = data.setdefault(w[0],[dict(),w[3]])[0].setdefault(w[1],dict()).setdefault(w[2],[list(),list()])
        vol_entry[0].append((w[7],w[8],w[9]))
        if len(data[w[0]]) == 2:
            data[w[0]].extend(w[4:7])
    return Dataset(data)

# Reads a CSV produced by write_csv and returns the contents as a Dataset object
def read_csv(in_file, seperator = ';', sub_seperator = '|'):
    data = dict()
    for l in in_file.readlines():
        if seperator == ',':    # again special cased
            w = re.findall('"([^"]*)"', l)
        else:
            w = l.replace('\n','').replace('"','').split(seperator)
        if len(w) == 10:    # pared manifolds weren't supported when this csv was written out
            w.append('')    # acceptable substitute
        vol_entry = data.setdefault(w[1],[dict(),w[5]])[0].setdefault(w[2],dict()).setdefault(w[4],[list(),list()])
        vol_entry[0].append((w[0],w[9],w[6]))
        vol_entry[1].extend(w[10].split(sub_seperator))
        vol_entry[1] = list(set(vol_entry[1]))  # remove duplicates
        if len(data[w[1]]) == 2:
            data[w[1]].extend([w[3],w[7],w[8]])
    return Dataset(data)

# Returns the list as a string with the given seperator and no brackets
def list_str(lst,sep):
    ret = ''
    for x in lst:
        ret += str(x)+sep
    return ret[:-1*len(sep)]    # remove extra seperator

# Writes a CSV file containing the mainfolds records as shown below.
# Note that pared manifolds are currently ignored.
def write_csv(out_file, dataset, seperator = ';', sub_seperator = '|', append=False):
    if not append:
        out_file.write('Name'+seperator+
                        'InvariantTraceField'+seperator+
                        'Root'+seperator+
                        'NumberOfComplexPlaces'+seperator+
                        'Volume'+seperator+
                        'InvariantTraceFieldDegree'+seperator+
                        'SolutionType'+seperator+
                        'Disc'+seperator+
                        'Factored'+seperator+
                        'Tetrahedra'+seperator+
                        'ParedManifolds'+'\n')
    for p in sorted(dataset.get_polys(), key=lambda poly: (int(dataset.get_degree(poly)), poly)):
        for r in dataset.get_roots(p):
            deg = dataset.get_degree(p)
            ncp = dataset.get_ncp(p)
            disc = dataset.get_disc(p)
            fact_disc = dataset.get_factored_disc(p)
            for v in dataset.get_volumes(p,r):
                for m in dataset.get_manifold_data(p,r,v):
                    out_file.write('"'+m[0]+'"'+seperator)
                    out_file.write('"'+p+'"'+seperator)
                    out_file.write('"'+r+'"'+seperator)
                    out_file.write('"'+ncp+'"'+seperator)
                    out_file.write('"'+v+'"'+seperator)
                    out_file.write('"'+deg+'"'+seperator)
                    out_file.write('"'+m[2]+'"'+seperator)
                    out_file.write('"'+disc+'"'+seperator)
                    out_file.write('"'+fact_disc+'"'+seperator)
                    out_file.write('"'+m[1]+seperator)
                    out_file.write('"'+list_str(dataset.get_pared_manifolds(p,r,v),sub_seperator).replace(' ','')+'"\n')

def quick_write_csv(dataset, filenm, seperator = ';', sub_seperator = '|', append = False):
    try:
        f = open(filenm,'w')
        write_csv(f, dataset, seperator = seperator, sub_seperator = sub_seperator, append = append)
    except:
        f.close()
        raise

#### For backwards compatability
def pare_all_volumes(data):
    """
    Deprecated.  Use data.pare_all_volumes() instead
    """
    data.pare_all_volumes()

def cull_all_volumes(data, epsilon = EPSILON):
    """
    Deprecated.  Use data.cull_all_volumes() instead
    """
    data.cull_all_volumes(epsilon)

def _span_guesses(data):
    spans = dict()
    for poly in data.get_polys():
        poly_dict = spans.setdefault(poly,dict())
        ncp = 0
        try:
            ncp = int(data.get_ncp(poly))
        except ValueError:
            print('Some kind of problem with ncp ' + str(data.get_ncp(poly)) + "\n")
        if ncp < 1:
            continue
        try:
            roots = set([_get_conjs(r) for r in data.get_roots(poly)])  # turn roots into equivalence classes of conjugation
            for root in roots:
                vols = list()
                for other in data.get_roots(poly):
                    if _up_to_conjugates(root,other):
                        vols.extend([(v,data.get_nice_manifold_name(poly,other,v)) for v in data.get_volumes(poly,other)])
                vols = [v for v in vols if gen.pari(str(v[0]) + ' > 0.9') ]
                if not vols:
                    continue
                try:
                    poly_dict[root] = find_span(vols, ncp)  # HERE
                except ValueError as ve:
                    poly_dict[root] = ("Error (" + str(ve) + ")", 0, "Error")
        except Exception as e:
            print(traceback.format_exc())
            print(str(e))
            pass
    return spans

def is_int(fl, epsilon = EPSILON):
    return fl % 1 < epsilon or 1 - (fl % 1) < epsilon

def quick_write_spans(in_filenames, out_filename, out_seperator = ';'):
    if type(in_filenames) is str: # support laziness
        in_filenames = [in_filenames]
    lines = []
    for f in in_filenames:
        fi = open(f, 'r')
        fi.readline()
        lines.extend([l for l in fi.readlines()])
        fi.close()
    d = read_raw_csv(lines)
    d.remove_non_geometric_elements()
    write_spans(out_filename, d, seperator = out_seperator)

def read_spans(fname, seperator = ';'):
    f = open(fname,'r')
    f.readline()
    spans = dict()
    for l in f.readlines():
        w = l.replace('"','').replace(' ','').strip('\n').split(seperator)  # whitespace can cause weird problems
        for i in [4,5,7,8,9,10]:
            try:
                w[i] = w[i][2:-2].split("','")    # convert string back to list of strings
            except IndexError:
                break
        spans.setdefault(w[0],dict())[w[2]] = w[4:]
    return spans

def write_spans(fname, dataset, seperator = ';'):
    d = dataset # lazy
    s = _span_guesses(d)
    f = open(fname, 'w')
    f.write('Polynomial' + seperator + 'NumberOfComplexPlaces' + seperator + 'Root' + seperator + 'SpanDimension' + seperator + 'VolumeSpan' + seperator + 'ManifoldSpan' + seperator + 'FitRatio\n')
    for p,pd in s.items():
        for r,re in pd.items():
            if str(re[1]) != '0':
                f.write('"' + str(p) + '"' + seperator)
                f.write('"' + str(d.get_ncp(p)) + '"' + seperator)
                f.write('"' + str(r) + '"' + seperator)
                f.write('"' + str(len(re[0])) + '"' + seperator)
                f.write('"' + str(re[0]) + '"' + seperator)
                f.write('"' + str(re[2]) + '"' + seperator)
                f.write('"' + str(re[1]) + '"\n')
    f.close()

# Manages interfacing with a dictionary containing the spans.
# There are two possible forms:
# poly : root : [[spanning_vols], fit_ratio, [spanning_names]]
# poly : root : [[spanning_vols], fit_ratio, [spanning_names], [good_pseudo(vols,names)], pseudo_fit_ratio, [bad_pseudo(vols,names)]]
# The latter form is used after deciding to fit some pseudovols (as a VolumeData) against a SpanData;
# Doing so produces a VolumeData of those pseudovols we just couldn't fit, which can be written out as usual
# Also, this seems prone (for some reason) to causing stack overflows in PARI
class SpanData:

    def __init__(self, data_dict, fails_dict = None):
        self.data = data_dict
        self.nice_fits = dict()
        if data_dict:
            p = data_dict.keys().__iter__().next()
            r = data_dict[p].keys().__iter__().next()  # this really shouldn't fail
            s = len(data_dict[p][r])
            if s != 3 and s != 6:
                raise ValueError    # input looks wack
            self.fitted = (s == 6)  # records if we are in the second form described above
        else:
            self.fitted = False
        self.fit_fails = None
        if fails_dict:
            self.fit_fails = fails_dict
        elif self.fitted:
            self.fit_fails = dict()
        for p in data_dict.keys():  # got to reformat
            for r in data_dict[p].keys():
                data_dict[p][r] = list(data_dict[p][r])

    def get_polys(self):
        return self.data.keys()

    def get_roots(self, poly):
        return self.data[poly].keys()

    def get_spans(self, poly, root):
        return self.data[poly][root][:3]

    def get_pseudo_data(self, poly, root):
        return self.data[poly][root][3:]

    # If fitted, gives a dict manifold --> itf --> root --> volume --> lindep
    # where lindep is the result of lindepping get_spans(ift,root) + [volume]
    def get_nice_fits(self):
        return self.nice_fits

    # Given a valid filename or open, write-ready file object, writes our data out to it.
    def write_to_csv(self, outfile, dset, seperator = ';', append = False):
        if type(outfile) == str:
            if append:
                f = open(outfile,'a')
            else:
                f = open(outfile,'w')
        else:
            f = outfile
        try:
            if not append:
                f.write('Polynomial' + seperator + 'NumberOfComplexPlaces' + seperator + 'Root' + seperator + 'SpanDimension' + seperator + 'VolumeSpan' + seperator + 'ManifoldSpan' + seperator + 'FitRatio')
                if self.fitted:
                    f.write(seperator + 'SolvedPseudoVolumes' + seperator + 'SolvedNames' + seperator + 'UnsolvedPseudoVolumes' + seperator + 'UnsolvedNames' + seperator + 'PseudoFitRatio')
                f.write('\n')
            for p in self.get_polys():
                for r in self.get_roots(p):
                    re = self.data[p][r]
                    f.write('"' + str(p) + '"' + seperator)
                    try:
                        f.write('"' + str(dset.get_ncp(p)) + '"' + seperator)
                    except: # don't give up because dset was surprised
                        f.write('"?"' + seperator)
                    f.write('"' + str(r) + '"' + seperator)
                    f.write('"' + str(len(re[0])) + '"' + seperator)
                    f.write('"' + str(re[0]) + '"' + seperator)
                    f.write('"' + str(re[2]) + '"' + seperator)
                    f.write('"' + str(re[1]) + '"')
                    if self.fitted:
                        f.write(seperator)
                        if len(re) == 6:
                            f.write('"' + str([t[0] for t in re[3]]) + '"' + seperator)
                            f.write('"' + str([t[1] for t in re[3]]) + '"' + seperator)
                            f.write('"' + str([t[0] for t in re[5]]) + '"' + seperator)
                            f.write('"' + str([t[1] for t in re[5]]) + '"' + seperator)
                            f.write('"' + str(re[4]) + '"')
                        else:
                            f.write('"None"' + seperator)
                            f.write('"None"' + seperator)
                            f.write('"None"' + seperator)
                            f.write('"None"' + seperator)
                            f.write('"1"')
                    f.write('\n')
        finally:
            if type(outfile) == str:
                f.close()

    def fit(self, voldata, maxcoeff = MAX_COEFF, max_ldp_tries = 4, max_itf_degree = MAX_ITF):
        def _fresz(p,r):    # if not already done, change data[p][r] to bigger format
            if len(self.data[p][r]) == 3:
                self.data[p][r].extend([list(),0,list()])
        def _fit(p,rec):      # this exists to break multiple layers
            cand = None     # previous best fit
            for tf in get_potential_trace_fields(p):
                tf = tf.replace(' ','')
                if tf in self.get_polys():
                    for r in self.get_roots(tf):
                        if 'Error' in self.data[tf][r]:
                            print 'Couldn\'t handle '+str(self.data[tf][r]) # can't handle the format
                            continue                                        # so we skip this one
                        ldp = _pari_lindep(self.get_spans(tf,r)[0]+[rec[0]], maxcoeff = maxcoeff, max_tries = max_ldp_tries)
                        if ldp and ldp[-1] != 0:
                            if abs(ldp[-1]) == 1:    # the match was perfect, update the data
                                _fresz(tf,r)
                                self.data[tf][r][3].append(rec)
                                self.nice_fits.setdefault(rec[1],dict()).setdefault(tf,dict()).setdefault(r,dict())[rec[0]] = ldp
                                return
                            else:   # the match was imperfect, maybe a better fit awaits
                                if not cand or cand[1][-1] > ldp[-1]:   # better than previous best fit
                                    cand = ((tf,r),ldp) # we store the whole lindep for later
            if cand:    # have a rational but not integral fit
                _fresz(cand[0][0],cand[0][1])
                self.data[cand[0][0]][cand[0][1]][5].append(rec)
                self.nice_fits.setdefault(rec[1],dict()).setdefault(cand[0][0],dict()).setdefault(cand[0][1],dict())[rec[0]] = cand[1]
            else:       # no rational fit, store the failure
                self.fit_fails.setdefault(p,list()).append(rec)
                # TODO store a trivial nicefits entry here
        if not self.fitted:
            self.fitted = True
            if not self.fit_fails:
                self.fit_fails = dict()
        for p in voldata.get_polys():
            for v in voldata.get_volumes(p):
                for m in voldata.get_manifolds(p,v):
                    _fit(p,(v,m))   # TODO fix this innefficent implementation (much redundnacy; should be once / v)
        for p in self.get_polys():      # got to recalc psuedo fit ratios
            for r in self.get_roots(p):
                if len(self.data[p][r]) == 6:    # only operate if pseudo volumes in play
                    if self.data[p][r][5]:  # TODO: make this record if coeff > fit ratio
                        dim = len(self.data[p][r][0])
                        vecs = list()
                        for n in xrange(dim):   # put in basis unit vectors
                            vecs.append([0]*dim)
                            vecs[n][n] = 1
                        for v in [rec[0] for rec in self.data[p][r][5]]:    # put in vectors for non-integral combinations
                            pldp = _pari_lindep(self.data[p][r][0]+[v])
                            vecs.append([Fraction(numerator = -1*x, denominator = pldp[-1]) for x in pldp[:-1]])
                        dets = list()
                        for c in combinations(vecs,dim):
                            dets.append(abs(det(c)))
                        self.data[p][r][4] = _gcd([Fraction(d) for d in dets])
                    else:   # all volumes fit integrally, so
                        self.data[p][r][4] = 1

    # Returns a dict poly : (volume, manifold) of manifolds that couldn't be fitted in the spans.
    def get_fit_failures(self):
        return self.fit_fails

    def write_failures(self, outfile, seperator = ';', append = False):
        if type(outfile) == str:
            if append:
                f = open(outfile,'a')
            else:
                f = open(outfile,'w')
        else:
            f = outfile
        try:
            if not append:
                f.write('TraceField' + seperator + 'Volume' + seperator + 'Manifold\n')
            for p in self.fit_fails.keys(): # was this working without keys?
                for rec in self.fit_fails[p]:
                    f.write('"'+str(p)+'"'+seperator)
                    f.write('"'+str(rec[0])+'"'+seperator)
                    f.write('"'+str(rec[1])+'"\n')
        finally:
            if type(outfile) == str:
                f.close()

    # Writes out the linear combinations producing exotic volumes in a relatively readable format as described below
    # The format for the combination is:
    # k*exotic_man=k*man_1+-k*man_2+-k*man_3...
    # where k are some nonzero integers (so no if one is negative), +- is + or -, exotic_man is Manifold,
    # and the other manifolds names stand in for their geometric volumes
    def write_nice_fits(self, outfile, seperator = ';', append = False):
        if type(outfile) == str:
            if append:
                f = open(outfile,'a')
            else:
                f = open(outfile,'w')
        else:
            f = outfile
        try:
            if not append:
                f.write('Manifold'+seperator+'InvTraceField'+seperator+'Root'+seperator+'Volume'+seperator+'Combination\n')
            for m in self.nice_fits.keys():
                for itf in self.nice_fits[m].keys():
                    for r in self.nice_fits[m][itf].keys():
                        for v in self.nice_fits[m][itf][r].keys():
                            ldp = self.nice_fits[m][itf][r][v]
                            comb = str(ldp[-1])+'*'+m+'='
                            for n in xrange(len(ldp)-1):
                                if n != 0 and -1*ldp[n] > 0:  # don't add a plus sign for the first term
                                    comb += '+'
                                if ldp[n] != 0:
                                    comb += str(-1*ldp[n])+'*'+self.get_spans(itf,r)[1][n]
                            f.write('"'+m+'"'+seperator)
                            f.write('"'+itf+'"'+seperator)
                            f.write('"'+r+'"'+seperator)
                            f.write('"'+v+'"'+seperator)
                            f.write('"'+comb+'"\n')
            for p in self.fit_fails.keys():
                for rec in self.fit_fails[p]:
                    f.write('"'+str(rec[1])+'"'+seperator)
                    f.write('"'+str(p)+'"'+seperator)
                    f.write('"'+'TraceField'+'"'+seperator)
                    f.write('"'+str(rec[0])+'"'+seperator)
                    f.write('"'+'None'+'"\n')
        finally:
            if type(outfile) == str:
                f.close()

# returns a SpanData for the given dataset
def get_data_object(dset):
    return SpanData(_span_guesses(dset))


# Accepts volumes as strings since we store them that way.
# Returns the dependancy found (if any) as a list of integers if all coefficents are <= maxcoeff or maxcoeff is nonpositive;
# otherwise, it returns []
def _pari_lindep(str_vols, maxcoeff = MAX_COEFF, max_tries = 50):
    vols = list(str_vols)   # in case someone sent some other type collection
    vec = None

    num_tries = 0
    while num_tries < max_tries:
        vec = str(pari(str(vols).replace("\'",'')).lindep())[1:-2].replace(' ','').split(',')
        num_tries += 1

    if not vec or vec == ['']: # no input
        #######
        print('Sometimes went wrong calculating ' + str(str_vols))
        #######

        return list()
    o = [int(v) for v in vec]
    if maxcoeff > 0:
        for x in o:
            if abs(x) > maxcoeff:
                o = list()
                break
    return o

def _lcm(a,b):
    return (a*b)/fractions.gcd(a,b) # gcd * lcm = a * b

# Hilariously, python's fractions.gcd does not accept Fractions or handle their floats correctly
def _gcd(fracts):
    denom = 1
    for q in fracts:
        denom = _lcm(denom, q.denominator)   # our denominator should be lcm of input denominators
    nums = [(q.numerator*denom)/q.denominator for q in fracts]  # create common denominator
    num = nums[0]
    for n in nums [1:]:
        nums = fractions.gcd(num,n) # our numerator should be gcd of input (same denom) numerators
    return Fraction(numerator = num, denominator = denom)

# Test code
if __name__ == '__main__':
    f = open('output.csv','r')
    f.readline() # skip header
    d = read_raw_csv_from_file(f)
    f.close()
    pare_all_volumes(d)
    cull_all_volumes(d)
    g = open('newoutput.csv','w')
    write_csv(g,d)
    g.close()
