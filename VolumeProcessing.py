#!/usr/bin/python

import re
import sys
import traceback
import fractions

from numpy.linalg import det
from SpanFinder import find_span
from PseudoVols import *
from cypari import *
from fractions import Fraction
from itertools import combinations

EPSILON = .0000000000001
MAX_COEFF = 4096

# This class is just a wrapper for the structure storing polynomial/volume data.
# Having it avoids opaque references to the particular way data is stored that might change in the future.

class dataset:
    def __init__(self, data_dict = dict()):
        self.data = data_dict

    def get_polys(self):
        return self.data.keys()

    def get_roots(self,poly):
        return self.data[poly][0].keys()

    def get_degree(self,poly):
        return self.data[poly][1]

    def get_ncp(self,poly):
        return self.data[poly][2]

    def get_disc(self,poly):
        return self.data[poly][3]

    def get_factored_disc(self,poly):
        return self.data[poly][4]

    def get_volumes(self,poly,root):
        return self.data[poly][0][root].keys()

    # checks if a given polynomial has a record in this dataset
    def has_poly(self,poly):
        return poly in self.data.keys()

    # returns a geometric manifold's record, or None failing that
    def get_geom_manifold(self,poly,root,vol):
        try:
            return self.data[poly][0][root][vol][0][0]
        except:
            return (None,None,None)

    # Returns triplets of a manifold's name, number of simplices, and solution type
    def get_manifold_data(self,poly,root,vol):
        return self.data[poly][0][root][vol][0]

    # Returns names of manifolds with vol that weren't chosen as representative
    def get_pared_manifolds(self,poly,root,vol):
        return self.data[poly][0][root][vol][1]

    # Removes a volume, returning its manifold data
    def remove_volume(self,poly,root,vol):
        rec = self.data[poly][0][root].get(vol)
        del self.data[poly][0][root][vol]
        return rec

    # Returns a dataset with the data from self and other; in case of a conflict, other's values beat self's or both are kept
    # Therefore, one is advised to use this on disjoint datasets or pare volumes afterwards
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
        return dataset(data_dict = new_data)

    # Return a triplet containing data for the manifold of smallest volume with the given field
    def get_representative_element(self, poly, root):
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
        newdata = dict()
        for poly in self.get_polys():
            newdata[poly] = [dict()]+self.data[poly][1:]    # initialize list of volumes to be empty
            for root in self.get_roots(poly):
                md = self.get_representative_element(poly,root)
                if md is not None:  # we actually have something geometric for this root
                    newdata[poly][0][root] = {md[0] : [[md[1]],list()]}
            if newdata[poly][0] == dict():  # no roots gave us geometric solutions
                del newdata[poly]
        return dataset(data_dict = newdata)

    # Returns false if contents look very wrong (no x in polynomial slot, etc.)
    # Only checks very shallowly for the first record data.__iter__.next() returns, so no guarantees
    def sane(self):
        try:
            try:
                p = self.data.keys().__iter__().next()
            except StopIteration:
                return True # empty dataset is sane
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

    # removes all manifolds from all volumes if they are anything but 'geometric'
    def remove_non_geometric_elements(self):
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

    # Currently returns the shortest choice; could be made more sophisticated
    def get_nice_manifold_name(self,poly,root,vol):
        nms = [rec[0] for rec in self.get_manifold_data(poly,root,vol)]
        nms.extend(self.get_pared_manifolds(poly,root,vol))
        opt = ('', sys.float_info.max)
        for nm in nms:
            if len(nm) < opt[1]:    # TODO: finish _niceness and replace len with it
                opt = (nm,len(nm))
        return opt[0]

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

    # If some volumes differ by less than epsilon, combine them, keeping one of them arbitrarily.
    # Slightly nondeterministic; may not get large, dense (compared to epsilon) clumps if iteration order is unfavorable.
    # But that is very unlikely to happen,
    # and when it is (lattice generators very close to each other) you don't want smush to work anyway
    # Only usefull if you don't want to cull;
    # Only briefly tested
    def smush_volumes(self, epsilon = EPSILON):
        balls = list()        
        for p in self.get_polys():
            for r in self.get_roots(p):
                vol_data = self.data[p][0][r]
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
                    

def quick_read_csv(filenm, seperator = ';', sub_seperator = '|'):
    try:    
        f = open(filenm,'r')
        f.readline()
        d = read_csv(f, seperator = seperator, sub_seperator = sub_seperator)
        f.close()
        return d
    except:
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
        w = [w[0],w[9],w[4],w[1],w[5],w[2],w[6],w[3],w[7],w[8]]
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
    return dataset(data)

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
    return dataset(data)

# Reads a CSV produced by write_csv and returns the contents as a dataset object
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
    return dataset(data)

# Returns the list as a string with the given seperator and no brackets
def list_str(lst,sep):
    ret = '['
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

# Removes redundant manifold records with the same trace field (and root) and volume
def pare_all_volumes(data):
    for p in data.get_polys():
        for r in data.get_roots(p):
            for v in data.get_volumes(p,r):
                pare_volume(data,p,r,v)

def pare_volume(data,poly,root,vol): # TODO: move these 4 methods into the class
    mdata = data.get_manifold_data(poly,root,vol)
    mpared = data.get_pared_manifolds(poly,root,vol)
    while len(mdata) > 1:
        mpared.append(mdata.pop(1)[0])

# Removes volumes that are integer multiples of another volume
def cull_all_volumes(data):
    for p in data.get_polys():
        for r in data.get_roots(p):
            cull_volumes(data,p,r)

def cull_volumes(data,poly,root):
    vols = data.get_volumes(poly,root)
    # vols = data.data[poly][0][root].keys()
    i = 0
    while i < len(vols) - 1:
        j = i + 1
        while j < len(vols):
            try:
                if is_int(float(vols[i])/float(vols[j])) and gen.pari(vols[i] + ' > ' + vols[j]):
                    # We have to throw away (culled) manifold names to let all culled manifolds have the same volume
                    # [j] divides [i] so remove [i]
                    data.remove_volume(poly,root,vols.pop(i))
                    # i is already effectivley incremented, so we must offset it
                    i = i-1
                    break
                elif is_int(float(vols[j])/float(vols[i])):
                    # this time, remove [j]
                    data.remove_volume(poly,root,vols.pop(j))
                    # j is effectivley incremented, no need to do it
                else:
                    j += 1
            except (ValueError, ZeroDivisionError): # bad quotient; not a linear combination either way so...
                j += 1
        i += 1

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
            for root in data.get_roots(poly):
                vols = [(gen.pari(v),data.get_geom_manifold(poly,root,v)[0]) for v in data.get_volumes(poly, root)]
                vols = [v for v in vols if gen.pari(str(v[0]) + ' > 0.9') ]
                if not vols:
                    continue
                try:
                    poly_dict[root] = find_span(vols, ncp)
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
    pare_all_volumes(d)
    cull_all_volumes(d)
    write_spans(out_filename, d, seperator = out_seperator)
    

def read_spans(fname, seperator = ';'):
    f = open(fname,'r')
    f.readline()
    spans = dict()
    for l in f.readlines():
        w = l.replace('"','').strip('\n').split(seperator)
        spans.setdefault(w[0],dict())[w[2]] = w[4:]    # with any luck, this is only set once
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
class SpanData:

    def __init__(self, data_dict, fails_dict = None):
        self.data = data_dict
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
                    f.write('"' + str(dset.get_ncp(p)) + '"' + seperator)
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

    def fit(self, voldata, maxcoeff = MAX_COEFF):
        def _fresz(p,r):    # if not already done, change data[p][r] to bigger format
            if len(self.data[p][r]) == 3:
                self.data[p][r].extend([list(),0,list()])
        def _fit(p,rec):      # this exists to break multiple layers
            cand = None     # previous best fit
            for tf in get_potential_trace_fields(p):
                tf = tf.replace(' ','')
                if tf in self.get_polys():
                    for r in self.get_roots(tf):
                        if 'Error' in self.data[tf][r]: # can't handle the format (TODO?)
                            continue                    # so we skip this one
                        ldp = _pari_lindep(self.get_spans(tf,r)[0]+[rec[0]], maxcoeff = maxcoeff)
                        if ldp and ldp[-1] != 0:
                            if abs(ldp[-1]) == 1:    # the match was perfect, update the data
                                _fresz(tf,r)
                                self.data[tf][r][3].append(rec)
                                return
                            else:   # the match was imperfect, maybe a better fit awaits
                                if not cand or cand[1] > ldp[-1]:   # better than previous best fit
                                    cand = ((tf,r),ldp[-1])
            if cand:    # have a rational but not integral fit
                _fresz(cand[0][0],cand[0][1])
                self.data[cand[0][0]][cand[0][1]][5].append(rec)
            else:       # no rational fit, store the failure
                self.fit_fails.setdefault(p,list()).append(rec)
        if not self.fitted:
            self.fitted = True
            if not self.fit_fails:
                self.fit_fails = dict()
        for p in voldata.get_polys():
            data = voldata.get_volume_data(p)
            for rec in data:
                _fit(p,rec)
        for p in self.get_polys():      # got to recalc psuedo fit ratios
            for r in self.get_roots(p):
                if len(self.data[p][r]) == 6:    # only operate if pseudo volumes in play
                    if self.data[p][r][5]:  # take gcd using basis vectors and pseudo vols; hopefully equivalent by linear algebra TODO: verify
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
            for p in self.fit_fails:
                for rec in self.fit_fails[p]:
                    f.write('"'+str(p)+'"'+seperator)
                    f.write('"'+str(rec[0])+'"'+seperator)
                    f.write('"'+str(rec[1])+'"\n')
        finally:
            if type(outfile) == str:
                f.close()

# returns a SpanData for the given dataset                                    
def get_data_object(dset):
    return SpanData(_span_guesses(dset))

# Accepts volumes as strings since we store them that way.
# Returns the dependancy found (if any) as a list of integers if all coefficents are <= maxcoeff or maxcoeff is nonpositive;
# otherwise, it returns []
def _pari_lindep(str_vols, maxcoeff = MAX_COEFF):
    vols = list(str_vols)   # in case someone sent some other type collection
    vec = str(pari(str(vols).replace('\'','')).lindep())[1:-2].replace(' ','').split(',')
    if not vec or vec == ['']: # no input; TODO implement more elegantly
        return list()
    o = [int(v) for v in vec]
    if maxcoeff > 0:
        for x in o:
            if x > maxcoeff:
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
