import copy
import re
import sys
import traceback
import fractions

from SpanFinder import find_span, find_borel_matrix
from PseudoVols import VolumeData, get_potential_trace_fields
from VolumeUtilities import *
from cypari import *
from fractions import Fraction
from itertools import combinations
from snappy import *

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
                    self.pare_volume(p,r,v)

    def pare_volume(self,poly,root,vol):
        """
        Compress a (polynomial, root, volume) triple, so that only one
        representative manifold is stored for it.  Manifolds which are
        discarded have their names stored, and can be retrieved by
        get_pared_volumes()
        """
        mdata = self.get_manifold_data(poly,root,vol)
        mpared = self.get_pared_manifolds(poly,root,vol)
        while len(mdata) > 1:
            mpared.append(mdata.pop(1)[0])

    def cull_all_volumes(self, epsilon = EPSILON):
        """
        Remove all volumes that are integer multiples of another,
        smaller volume, where each are for the same polynomial and
        root. These volumes are not pared (and so cannot be retrieved by
        get_pared_volumes()), they are removed outright from the
        Dataset. For large Datasests, this may free resources and make
        dealing with the Dataset faster.
        """
        for p in self.get_polys():
            for r in self.get_roots(p):
                self.cull_volumes(p,r,epsilon = epsilon)

    def cull_volumes(self,poly,root,epsilon = EPSILON):
        """
        Remove all volumes that are integer multiples of another,
        smaller volume, where each satisfy the given polynomial and
        root. These volumes are not pared (and so cannot be retrieved by
        get_pared_volumes()), they are removed outright from the
        Dataset. For large Datasests, this may free resources and make
        dealing with the Dataset faster.
        """
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

    def combine_with(self,other):
        """
        Returns a Dataset with the merged contents of other; in case of
        a conflict, other's values take priority, though both are kept
        if there would be no conflict (as should be the case for most
        Datasets).

        Therefore, it is advised to use this on disjoint Datasets or
        pare volumes afterwards.
        """
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
        volume with the given field.
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

    def quick_write_csv(self, filenm, separator = ';', sub_separator = '|', append = False):
        """
        Write out the Dataset to an output file.  If append is set to
        False, the file, if extant, is overwritten, if not, the file is
        assumed to be complete and well-formed, including a header.
        """
        try:
            f = open(filenm,'w')
            write_csv(f, dataset, separator = separator, sub_separator = sub_separator, append = append)
        except:
            f.close()
            raise


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

    def smush_volumes(self, epsilon = EPSILON):
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
        """
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

    def write_spans(self, fname, separator = ';', skip_borel = False, borel_shape_field_degree = 32):
        """
        Collect the manifolds of this dataset into spanning lattices,
        writing the results out to the file specified by fname.

        Note: Computing borel regulators is quite intensive. If time
        is at all a concern and Borel regulators are not all desired,
        setting skip_borel to True will speed this up by many orders
        of magnitude. The parameter borel_shape_field_degree is passed
        directly to snap as by the command `set degree'. As it
        increases, the computations become slower, but more results
        are obtained.
        """
        s = _span_guesses(self)
        f = open(fname, 'w')
        f.write('Polynomial' + separator +
                'Degree' + separator +
                'NumberOfComplexPlaces' + separator +
                'Root' + separator +
                'SpanDimension' + separator +
                'VolumeSpan' + separator +
                'ManifoldSpan' + separator +
                'FitRatio' + separator +
                'BorelRegulatorMatrix' + separator +
                'BorelRegulatorDeterminant\n')
        for p,pd in s.items():
            for r,re in pd.items():
                if str(re[1]) != '0':
                    borel_regs = 'N/A'
                    borel_det = None
                    if skip_borel:
                        borel_regs = 'Not computed'
                    else:
                        try:
                            borel_regs, borel_det = find_borel_matrix(re[2],
                                                                      borel_shape_field_degree)
                        except:
                            pass

                    if not borel_det:
                        borel_det = 'N/A'

                    f.write('"' + str(p) + '"' + separator)
                    f.write('"' + str(gen.pari(p).poldegree()) + '"' + separator)
                    f.write('"' + str(self.get_ncp(p)) + '"' + separator)
                    f.write('"' + str(r) + '"' + separator)
                    f.write('"' + str(len(re[0])) + '"' + separator)
                    f.write('"' + str(re[0]) + '"' + separator)
                    f.write('"' + str(re[2]) + '"' + separator)
                    f.write('"' + str(re[1]) + '"' + separator)
                    f.write('"' + str(borel_regs) + '"' + separator)
                    f.write('"' + str(borel_det) + '"\n')
        f.close()

    def search_for_manifold(self,man):
        """
        Given a valid Manifold object (or a valid Manifold name, return
        information on it.

        Returned data is of the form [InvariantTraceField, Root, Volume,
        SolutionType, GeomAlternative, NiceAlternative]

        Where all elements are strings.  Note that GeomAlternative and
        NiceAlternative may very well be the same manifold.
        """
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

class SpanData:
    """
    Manages interfacing with a dictionary containing the spans.  There are
    two possible forms:

    poly : root : [[spanning_vols], fit_ratio, [spanning_names]]

    poly : root : [[spanning_vols], fit_ratio, [spanning_names],
    [good_pseudo(vols,names)], pseudo_fit_ratio, [bad_pseudo(vols,names)]]

    The latter form is used after deciding to fit some pseudovols (as a
    VolumeData) against a SpanData; Doing so produces a VolumeData of those
    pseudovols we just couldn't fit, which can be written out as usual.

    Also, this seems prone (for some reason) to causing stack overflows
    in PARI
    """
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
        """
        Return a list of all polynomials represeted in this SpanData, as
        strings
        """
        return self.data.keys()

    def get_roots(self, poly):
        """
        Return a list of all roots of the given polynomial (input as a
        string) in this SpanData, as strings
        """
        return self.data[poly].keys()

    def get_spans(self, poly, root):
        """
        Return the spans for the given polynomial and root (input as
        strings) of this SpanData
        """
        return self.data[poly][root][:3]

    def get_pseudo_data(self, poly, root):
        """
        Return the pseudovolume data for the given polynomial and root
        (input as strings) of this SpanData
        """
        return self.data[poly][root][3:]

    def get_nice_fits(self):
        """
        If fitted, gives a dictionary chain.  The key progression of this chain is

        manifold name
          invariant trace field
            root
              volume
                the result of pari's lindep() on the span for this
                polynomial and root and this volume.

        In the lindep result, the volume was inserted into the linear
        dependence as the last element.
        """
        return self.nice_fits

    def write_to_csv(self, outfile, dset, separator = ';', append = False):
        """
        Write these span results out to outfile as a csv.
        """
        if type(outfile) == str:
            if append:
                f = open(outfile,'a')
            else:
                f = open(outfile,'w')
        else:
            f = outfile
        try:
            if not append:
                f.write('Polynomial' + separator + 'NumberOfComplexPlaces' + separator + 'Root' + separator + 'SpanDimension' + separator + 'VolumeSpan' + separator + 'ManifoldSpan' + separator + 'FitRatio')
                if self.fitted:
                    f.write(separator + 'SolvedPseudoVolumes' + separator + 'SolvedNames' + separator + 'UnsolvedPseudoVolumes' + separator + 'UnsolvedNames' + separator + 'PseudoFitRatio')
                f.write('\n')
            for p in self.get_polys():
                for r in self.get_roots(p):
                    re = self.data[p][r]
                    f.write('"' + str(p) + '"' + separator)
                    try:
                        f.write('"' + str(dset.get_ncp(p)) + '"' + separator)
                    except: # don't give up because dset was surprised
                        f.write('"?"' + separator)
                    f.write('"' + str(r) + '"' + separator)
                    f.write('"' + str(len(re[0])) + '"' + separator)
                    f.write('"' + str(re[0]) + '"' + separator)
                    f.write('"' + str(re[2]) + '"' + separator)
                    f.write('"' + str(re[1]) + '"')
                    if self.fitted:
                        f.write(separator)
                        if len(re) == 6:
                            f.write('"' + str([t[0] for t in re[3]]) + '"' + separator)
                            f.write('"' + str([t[1] for t in re[3]]) + '"' + separator)
                            f.write('"' + str([t[0] for t in re[5]]) + '"' + separator)
                            f.write('"' + str([t[1] for t in re[5]]) + '"' + separator)
                            f.write('"' + str(re[4]) + '"')
                        else:
                            f.write('"None"' + separator)
                            f.write('"None"' + separator)
                            f.write('"None"' + separator)
                            f.write('"None"' + separator)
                            f.write('"1"')
                    f.write('\n')
        finally:
            if type(outfile) == str:
                f.close()

    def fit(self, voldata, maxcoeff = MAX_COEFF, max_ldp_tries = MAX_LDP_TRIES, max_itf_degree = MAX_ITF):
        """
        Given a VolumeData object (from PseudoVols) representing some exotic volumes, this method attempts to see if we can generate them
        as linear combinations of the volumes in the spans. After calling this, you can write out the fits with write_nice_fits, and
        if you write out the SpanData object, data on the fits will be included.

        When this code is run, you will get a lot of "***   polynomial not in Z[X] in galoisinit." printed out; don't worry, that's normal.
        """
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
        """
        Returns a dictionary.  The keys in the dictionary are
        polynomials (as strings), the values are lists of tuples of
        (volume, manifold), for each volume and its accompanying
        manifold that have the invariant trace field of the given
        polynomial, but could not be fitted.
        """
        return self.fit_fails

    def write_failures(self, outfile, separator = ';', append = False):
        """
        Write the fit failures (see get_fit_failures) out to outfile as a csv
        """
        if type(outfile) == str:
            if append:
                f = open(outfile,'a')
            else:
                f = open(outfile,'w')
        else:
            f = outfile
        try:
            if not append:
                f.write('TraceField' + separator + 'Volume' + separator + 'Manifold\n')
            for p in self.fit_fails.keys(): # was this working without keys?
                for rec in self.fit_fails[p]:
                    f.write('"'+str(p)+'"'+separator)
                    f.write('"'+str(rec[0])+'"'+separator)
                    f.write('"'+str(rec[1])+'"\n')
        finally:
            if type(outfile) == str:
                f.close()


    def write_nice_fits(self, outfile, separator = ';', append = False):
        """
        Writes out the linear combinations producing exotic volumes in a
        relatively readable format as described below.

        The format for the combination is:


        k1 * exotic_man = k2 * man_1 +- k3 * man_2 +- k4 * man_3...

        where ki are each some nonzero integers (so no if one is
        negative), +- is + or -, exotic_man is Manifold, and the other
        manifolds names stand in for their geometric volumes.
        """
        if type(outfile) == str:
            if append:
                f = open(outfile,'a')
            else:
                f = open(outfile,'w')
        else:
            f = outfile
        try:
            if not append:
                f.write('Manifold'+separator+'InvTraceField'+separator+'Root'+separator+'Volume'+separator+'Combination\n')
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
                            f.write('"'+m+'"'+separator)
                            f.write('"'+itf+'"'+separator)
                            f.write('"'+r+'"'+separator)
                            f.write('"'+v+'"'+separator)
                            f.write('"'+comb+'"\n')
            for p in self.fit_fails.keys():
                for rec in self.fit_fails[p]:
                    f.write('"'+str(rec[1])+'"'+separator)
                    try:
                        f.write('"'+str(pari(str(p)).polredabs())+'"'+separator)
                    except: # cypari.gen.error or w/e from polredabs failing
                        f.write('"'+str(p)+'"'+separator)   # be consistent with get_potential_trace_field fail behaviour
                    f.write('"'+'TraceField'+'"'+separator)
                    f.write('"'+str(rec[0])+'"'+separator)
                    f.write('"'+'None'+'"\n')
        finally:
            if type(outfile) == str:
                f.close()

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

    # Dehn surgery penalties:
    in_dehn_fill = 0
    for l in nm[1:]:
        if l == '(':
            in_dehn_fill = 1
        elif l == ')':
            in_dehn_fill = 0
        elif l in '123456789' and in_dehn_fill:
            n += 0.5

    return n

def quick_read_csv(filenm, separator = ';', sub_separator = '|'):
    """
    Read in a csv (with header) and return a Dataset object representing
    it.  The csv should be in the form exected by read_csv
    """
    try:
        f = open(filenm,'r')
        f.readline()
        d = read_csv(f, separator = separator, sub_separator = sub_separator)
        f.close()
        return d
    except:
        if f:
            f.close()
        raise

# combines two output files from this program
def quick_combine_files(filenms, fileseps, out_filenm, out_separator = ';', out_append = False):
    """
    Given a list of filenames and file separators for each file, combine
    the volumes of each file into one output.  The result is similar to
    `cat`, in addition to handling the headers and conversions of
    separators, if desired.
    """
    dsets = list()
    for i in xrange(len(filenms)):
        inf = open(filenms[i],'r')
        try:
            inf.readline()  # skip header
            dsets.append(read_csv(inf, separator = fileseps[i]))
        finally:
            inf.close()
    for d in dsets[1:]:
        dsets[0].combine_with(d)
    if out_append:
        ouf = open(out_filenm, 'a')
    else:
        ouf = open(out_filenm, 'w')
    write_csv(ouf, dsets[0], separator = out_separator, append = out_append)

def quick_preprocess(in_filenm, out_filenm, in_separator = ';', out_separator = ';', out_append = False):
    """
    A convenience method to do the following: Read from the input
    filename, pare and cull the resulting dataset, then write it out ot
    the output filename.
    """
    inf = open(in_filenm,'r')
    try:
        inf.readline()  # skip header
        d = read_raw_csv(inf, separator = in_separator)
    finally:
        inf.close()
    pare_all_volumes(d)
    cull_all_volumes(d)
    if out_append:
        ouf = open(out_filenm,'a')
    else:
        ouf = open(out_filenm,'w')
    try:
        write_csv(ouf, d, separator = out_separator, append = out_append)
    finally:
        ouf.close()

def read_raw_csv_from_file(in_file, separator = ';'):
    """
    Read raw csv data in.  A header is not expected.  The file should be
    in the same format as expected by read_csv
    """
    return read_raw_csv(in_file.readlines(), separator)

def _up_to_conjugates(z,w):
    """
    Returns true if the given strings are equal or complex conjugates as
    formatted by snap: a+b*I, a-b*I
    """
    zp = re.findall(r'([+-]?[\d.]+)',z)
    wp = re.findall(r'([+-]?[\d.]+)',w)
    return len(zp) == len(wp) == 2 and zp[0] == wp[0] and up_to_sign(zp[1],wp[1])

def up_to_sign(x,y):
    """
    Returns true if one of the strings is just -the other. This method
    should only be applied to non sci-notation floats' strings
    """
    return re.search(r'[\d.]+',x).group() == re.search(r'[\d.]+',y).group()

def _get_conjs(z):
    """
    Given a+b*I (in string form), return a\pm b*I, where \pm is the
    unicode plus/minus character.
    """
    return z[:1]+z[1:].replace('+','\xb1').replace('-','\xb1')

def read_raw_csv(contents, separator = ';'):
    """
    Read in csv of the form

    Name;InvariantTraceField;Root;NumberOfComplexPlaces;Volume;InvariantTraceFieldDegree;SolutionType;Disc;DiscFactors;Tetrahedra
    """
    data = dict()
    # Obviously this code is highly sensative to any changes in the output format of VolumeFinder.py
    for l in contents:
        l = l.replace(' ', '')

        if separator == ',':    # special cased since ',' appears in Dehn surgery
            w = re.findall('"([^"]*)"', l)
        else:
            w = l.replace('\n','').replace('"','').split(separator)

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

def read_old_csv(in_file, separator = ';'):
    """
    Reads a CSV produced by write_csv and returns the contents as a
    dataset object.  This variant handles csvs before we swapped column
    order around a bit.
    """
    data = dict()
    for l in in_file.readlines():
        if separator == ',':    # again special cased
            w = re.findall('"([^"]*)"', l)
        else:
            w = l.replace('\n','').replace('"','').split(separator)
        vol_entry = data.setdefault(w[0],[dict(),w[3]])[0].setdefault(w[1],dict()).setdefault(w[2],[list(),list()])
        vol_entry[0].append((w[7],w[8],w[9]))
        if len(data[w[0]]) == 2:
            data[w[0]].extend(w[4:7])
    return Dataset(data)

def read_csv(in_file, separator = ';', sub_separator = '|'):
    """
    Read in a csv (without header) and return a Dataset object
    representing it.  The csv should be in the following form:

    Name;InvariantTraceField;Root;NumberOfComplexPlaces;Volume;InvariantTraceFieldDegree;SolutionType;Disc;DiscFactors;Tetrahedra
    """
    data = dict()
    for l in in_file.readlines():
        if separator == ',':    # again special cased
            w = re.findall('"([^"]*)"', l)
        else:
            w = l.replace('\n','').replace('"','').split(separator)
        if len(w) == 10:    # pared manifolds weren't supported when this csv was written out
            w.append('')    # acceptable substitute
        vol_entry = data.setdefault(w[1],[dict(),w[5]])[0].setdefault(w[2],dict()).setdefault(w[4],[list(),list()])
        vol_entry[0].append((w[0],w[9],w[6]))
        vol_entry[1].extend(w[10].split(sub_separator))
        vol_entry[1] = list(set(vol_entry[1]))  # remove duplicates
        if len(data[w[1]]) == 2:
            data[w[1]].extend([w[3],w[7],w[8]])
    return Dataset(data)

def _list_str(lst,sep):
    """
    Returns the list as a string with the given separator and no
    brackets.
    """
    ret = ''
    for x in lst:
        ret += str(x)+sep
    return ret[:-1*len(sep)]    # remove extra separator

def write_csv(out_file, dataset, separator = ';', sub_separator = '|', append=False):
    """
    Writes a CSV file containing the mainfolds records as shown
    below. Note that pared manifolds are currently ignored.
    """
    if not append:
        out_file.write('Name'+separator+
                        'InvariantTraceField'+separator+
                        'Root'+separator+
                        'NumberOfComplexPlaces'+separator+
                        'Volume'+separator+
                        'InvariantTraceFieldDegree'+separator+
                        'SolutionType'+separator+
                        'Disc'+separator+
                        'Factored'+separator+
                        'Tetrahedra'+separator+
                        'ParedManifolds'+'\n')
    for p in sorted(dataset.get_polys(), key=lambda poly: (int(dataset.get_degree(poly)), poly)):
        for r in dataset.get_roots(p):
            deg = dataset.get_degree(p)
            ncp = dataset.get_ncp(p)
            disc = dataset.get_disc(p)
            fact_disc = dataset.get_factored_disc(p)
            for v in dataset.get_volumes(p,r):
                for m in dataset.get_manifold_data(p,r,v):
                    out_file.write('"'+m[0]+'"'+separator)
                    out_file.write('"'+p+'"'+separator)
                    out_file.write('"'+r+'"'+separator)
                    out_file.write('"'+ncp+'"'+separator)
                    out_file.write('"'+v+'"'+separator)
                    out_file.write('"'+deg+'"'+separator)
                    out_file.write('"'+m[2]+'"'+separator)
                    out_file.write('"'+disc+'"'+separator)
                    out_file.write('"'+fact_disc+'"'+separator)
                    out_file.write('"'+m[1]+separator)
                    out_file.write('"'+_list_str(dataset.get_pared_manifolds(p,r,v),sub_separator).replace(' ','')+'"\n')

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

def quick_write_csv(data, filenm, separator = ';', sub_separator = '|', append = False):
    """
    Deprecated.  Use data.quick_write_csv() instead
    """
    data.quick_write_csv(filenm, separator, sub_separator, append)

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

def is_int(fl, epsilon = EPSILON): #TODO: move into utility methods
    """
    Return true if the floating point part of fl is within epsilon of 0.
    """
    return fl % 1 < epsilon or 1 - (fl % 1) < epsilon

def quick_write_spans(in_filenames, out_filename, out_separator = ';', skip_borel = False, borel_shape_field_degree = 32):
    """
    Compress input filenames

    Note: Computing borel regulators is quite intensive. If time is at
    all a concern and Borel regulators are not all desired, setting
    skip_borel to True will speed this up by many orders of
    magnitude. See write_spans for the optional
    borel_shape_field_degree argument.
    """
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
    d.write_spans(out_filename, separator = out_separator, skip_borel = skip_borel)

def read_spans(fname, separator = ';'):
    """
    Read in a span file, of the form
    Polynomial;NumberOfComplexPlaces;Root;SpanDimension;VolumeSpan;ManifoldSpan;FitRatio

    Returns a dictionary object (certainly NOT a Dataset) such that they
    keys are polynomials, and the values are dictionaries. These
    dictionaries have keys of roots and the values are [SpanDimension,
    VolumeSpan, ManifoldSpan, FitRatio.
    """
    f = open(fname,'r')
    f.readline()
    spans = dict()
    for l in f.readlines():
        w = l.replace('"','').replace(' ','').strip('\n').split(separator)  # whitespace can cause weird problems
        for i in [4,5,7,8,9,10]:
            try:
                w[i] = w[i][2:-2].split("','")    # convert string back to list of strings
            except IndexError:
                break
        spans.setdefault(w[0],dict())[w[2]] = w[4:]
    return spans


def get_spandata(dset):
    """
    Returns best guesses for spans generated by the given dataset
    """
    return SpanData(_span_guesses(dset))

def _pari_lindep(str_vols, maxcoeff = MAX_COEFF, max_tries = MAX_LINDEP_TRIES):
    """
    Given str_volumes, a list of volumes in string form, returns the
    dependancy found (if any) as a list of integers if all coefficents
    are <= maxcoeff or maxcoeff is nonpositive; otherwise, it returns []
    """
    vols = list(str_vols)   # in case someone sent some other type collection
    vec = None

    num_tries = 0
    while num_tries < max_tries:
        vec = str(pari(str(vols).replace("\'",'')).lindep(LINDEP_PRECISION))[1:-2].replace(' ','').split(',')
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
