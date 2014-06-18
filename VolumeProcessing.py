#!/usr/bin/python

import re

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

    # Returns pairs of a manifold's name and number of simplices
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

# Load a CSV file organized by manifold and reorganize it by polynomial and volume.
# The result: dict poly ---> (dict roots ----> (dict vols ---> (list (manifold name, tetrahedra), list (pared manifolds)), degree etc.)
def read_raw_csv(in_file):
    data = dict()
    # Obviously this code is highly sensative to any changes in the output format of VolumeFinder.py
    for l in in_file.readlines():
            # To avoid breaking this, make sure to quote all data fields.
            # Second replace is a bit of a hack b/c first was failing unexpectedly
            # w = l.replace('\n','').strip('"').replace('","','$').split('$')
            w = re.findall('"([^"]*)"', l) # Revert if you wish, I'm concerned about unknown, weird manifold names with '$' in them
            # Incase the disc was 1, a temporary hack:
            if len(w) == 8:
                w.append('')
            # w[0]: manifold name ---------------------------> m[0] for m in data[poly][0][root][vol][0]
            # w[1]: manifold simplices ----------------------> m[1] for m in data[poly][0][root][vol][0]
            # w[2]: volume ----------------------------------> v in data[poly][0][root].keys()
            # w[3]: invariant trace field polynomial --------> p in data.keys()
            # w[4]: polynomial degree -----------------------> data[poly][1]
            # w[5]: polynomial root -------------------------> r in data[poly][0].keys()
            # w[6]: polynomial number of complex places -----> data[poly][2]
            # w[7]: polynomial discriminant -----------------> data[poly][3]
            # w[8]: polynomial discriminant (factorized) ----> data[poly][4]
            # vr = data.setdefault(w[3],[dict(),w[4]])[0].setdefault(w[2],[list(),list(),w[5]])[0].append(w[0:2]) # OLD
            # # why was vr set just now and not used?
            vol_entry = data.setdefault(w[3],[dict(),w[4]])[0].setdefault(w[5],dict()).setdefault(w[2],[list(),list()])
            vol_entry[0].append((w[0],w[1]))
            if len(data[w[3]]) == 2:
                data[w[3]].extend(w[6:9])
            # print data[w[3]][1:] # DEBUG
    return dataset(data)

# Writes a CSV file containing the mainfolds records as shown below.
# Note that pared manifolds are currently ignored.
def write_csv(out_file, dataset, append=False):
    if not append:
        out_file.write('InvariantTraceField,Root,Volume,InvariantTraceFieldDegree,NumberOfComplexPlaces,Disc,Factored,Name,Tetrahedra\n')
    for p in dataset.get_polys():
        for r in dataset.get_roots(p):
            for v in dataset.get_volumes(p,r):
                for m in dataset.get_manifold_data(p,r,v):
                    out_file.write('"'+p+'",')
                    out_file.write('"'+r+'",')
                    out_file.write('"'+v+'",')
                    out_file.write('"'+dataset.get_degree(p)+'",')
                    out_file.write('"'+dataset.get_ncp(p)+'",')
                    out_file.write('"'+dataset.get_disc(p)+'",')
                    out_file.write('"'+dataset.get_factored_disc(p)+'",')
                    out_file.write('"'+m[0]+'",')
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
    # vols = data.get_volumes(poly,root)
    vols = data.data[poly][0][root].keys()
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
