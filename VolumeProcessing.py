#!/usr/bin/python

import re
import sys

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

    # Combines this dataset with the dataset other; in case of a difference, other's values beat self's
    def combine_with(self,other):
        self.update(other)

    # Return a triplet containing data for the manifold of smallest volume with the given field
    def get_representative_element(self, poly, root):
        minvol = (None, float(sys.maxint))       
        for v in self.get_volumes(poly,root):
            try:
                if float(v) < minvol[1] and 'geometric' in [rec[2] for rec in self.get_manifold_data(poly, root, v)]:
                    minvol = (v, float(v))
            except ValueError:
                continue    # probably v=0+; doesn't really matter which we pick anyway
        if minvol[0] is None:
            return None # means no geometric solutions were found for this field
        else:
            for m in self.get_manifold_data(poly,root,minvol[0]):
                if m[2] == 'geometric':
                    return (minvol[0],m)
                else:
                    print('Warning: this should never be reached.') # DEBUG
                    return None

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
def read_raw_csv(in_file, seperator = ';'):
    data = dict()
    # Obviously this code is highly sensative to any changes in the output format of VolumeFinder.py
    for l in in_file.readlines():
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
            # print data[w[3]][1:] # DEBUG
    return dataset(data)

# Reads csv from before the formatting change.
def read_old_raw_csv(in_file, seperator = ';'):
    data = dict()
    # Obviously this code is highly sensative to any changes in the output format of VolumeFinder.py
    for l in in_file.readlines():
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
            # print data[w[3]][1:] # DEBUG
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
def read_csv(in_file, seperator = ';'):
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


# Writes a CSV file containing the mainfolds records as shown below.
# Note that pared manifolds are currently ignored.
def write_csv(out_file, dataset, seperator = ';', append=False):
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
                        'Tetrahedra\n')
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
                    out_file.write('"'+m[1]+'"\n')

# Removes redundant manifold records with the same trace field (and root) and volume
def pare_all_volumes(data):
    for p in data.get_polys():
        for r in data.get_roots(p):
            for v in data.get_volumes(p,r):
                pare_volume(data,p,r,v)

def pare_volume(data,poly,root,vol):
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
                if is_int(float(vols[i])/float(vols[j])):
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

def is_int(fl, epsilon = .0000000000001):
    return fl % 1 < epsilon or 1 - (fl % 1) < epsilon

# Test code
if __name__ == '__main__':
    f = open('output.csv','r')
    f.readline() # skip header
    d = read_raw_csv(f)
    f.close()
    pare_all_volumes(d)
    cull_all_volumes(d)
    g = open('newoutput.csv','w')
    write_csv(g,d)
    g.close()
