#!/usr/bin/python

from snappy import *

# Returns a VolumeData object for manifolds with the given names
# Volumes' precision is based on pari, so set it there
def get_volume_data(man_nms, engine = 'magma', retrieve = True):
    recs = dict()
    for nm in man_nms:
        try:
            var = Manifold(nm).ptolemy_variety(2,'all')
            try:
                if retrieve:
                    sols = var.retrieve_solutions()
                else:
                    raise Exception('Coding too lazy!') 
            except:
                sols = var.compute_solutions(engine = engine)
            data = [(c.number_field(),c.volume_numerical()) for c in sols.flatten()]
            for d in data:
                for v in d[1]:
                    recs.setdefault(str(d[0]),list()).append((str(v),nm))
            for k in recs.keys():   # remove duplicates
                recs[k] = list(set(recs[k]))
        except Exception as e:
            print(str(e))+'; skipping '+nm
            continue
    return VolumeData(data = recs)

def get_potential_trace_fields(poly):
    pol = pari(poly)
    try:
        return [str(rec[0]) for rec in pol.nfsubfields()[1:] if _binmiss(rec[0].poldegree(),pol.poldegree())]   # poldegree returns int
    except: # we want cypari.gen.PariError, but no idea how to reference; fortunately, anything else will just raise again
        pol = pol.polredabs()
        return get_potential_trace_fields(str(pol))

def _binmiss(s,l):
    while s < l:
        s *= 2
    return s == l

# Wrapper for manipulating data on pseudo-volumes
class VolumeData:

    # structure: dict poly ---> list of (volume, manifold, ?fitted)
    def __init__(self, data = dict()):
        self.data = data

    def get_polys(self):
        return self.data.keys()

    def get_volumes(self,poly):
        return [rec[0] for rec in self.data[poly]]

    def get_manifolds(self,poly):
        return [rec[1] for rec in self.data[poly]]

    def get_volume_data(self,poly):
        return self.data[poly]

    # returns a VolumeData object containing the data from this and other; in case of a conflict, other's data takes precendence
    def combine_with(self,other):
        new_data = dict(self.data)
        for p in other.get_polys():
            new_data.setdefault(p,dict())
            new_data[p] = list(set(new_data[p].extend(other.get_volume_data(p))))
        return VolumeData(data = new_data)

    # given an (open, ready to write) file object or valid filename, writes the data
    def write_to_csv(self, output_file, seperator = ';', append = False):
        try:
            if type(output_file) == str:
                f = open(output_file,'w')
            else:
                f = output_file
            f.write('"TraceField"'+seperator+'"Volume"'+seperator+'"Manifold"'+'\n')
            for p in self.get_polys():
                for rec in self.get_volume_data(p):
                    f.write('"'+p+'"'+seperator)
                    f.write('"'+rec[0]+'"'+seperator)
                    f.write('"'+rec[1]+'"\n')
        finally:
            if type(output_file) == str and f:
                f.close()


    # given an (open, ready to read data) file object or valid filename, reads the file and returns a VolumeData that would write it
def read_volumedata_csv(infile, seperator = ';'):
    try:
        if type(infile) ==  str:
            f = open(infile,'r')
            f.readline()
        else:
            f = infile
        data = dict()
        for l in f.readlines():
            w = l.strip('\n').replace('"','').split(seperator)
            data.setdefault(w[0],list()).append(w[1:])
        return VolumeData(data = data)
    finally:
        if type(infile) == str and f:
            f.close()
