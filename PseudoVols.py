#!/usr/bin/python

from snappy import *
from ManifoldIterators import *

# Until this gets globalized
pari.set_real_precision(100)

# The same as calling get_volume_data(mans).write_to_csv(ofilenm) with the given parameters,
# except output will be written out every period manifolds and logs generated, instead of all at once.
def prepare_pvolume_file(maniter, ofilenm, append = False, engine = 'magma', retrieve = True, period = 100, seperator = ';'):
    ctr = 0
    block = list()
    done = False
    try:
        if append:
            f = open(ofilenm,'a')
        else:
            f = open(ofilenm,'w')
        while True:
            try:
                block.append(maniter.next())
                ctr += 1
            except StopIteration:
                done = True
            if ctr == period or done:
                print 'Processing '+str(block[0])+' to '+str(block[-1])+'.' 
                get_volume_data(ForwardingIterator(block.__iter__(),lambda m : str(m))).write_to_csv(f,append=append,seperator=seperator)
                append = True  # we must be appending after the first time
                ctr = 0
                block = list()
                if done:
                    break
    finally:
        f.close()

# Returns a VolumeData object containing exotic volumes for manifolds with the given names
# Volumes' precision is based on pari, so set it there
# set engine = None to skip computing the volumes 
def get_volume_data(man_nms, engine = 'magma', retrieve = True):
    recs = dict()
    for nm in man_nms:
        try:
            sols = None
            var = Manifold(nm).ptolemy_variety(2,'all')
            try:
                if retrieve:
                    sols = var.retrieve_solutions()
                else:
                    raise Exception('Coding too lazy!') 
            except Exception as e: # try using engine
                if engine:
                    sols = var.compute_solutions(engine = engine)
                else:
                    print(str(e))+'; skipping '+nm
            if sols:
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

    # structure: dict poly ---> dict volume ---> [manifolds]
    def __init__(self, data = dict()):
        self.data = data

    def get_polys(self):
        return self.data.keys()

    def get_volumes(self,poly):
        return self.data[poly].keys()

    def get_manifolds(self,poly,volume):
        return self.data[poly][volume]

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
                if append:
                    f = open(output_file,'a')
                else:
                    f = open(output_file,'w')
            else:
                f = output_file
            if not append:
                f.write('"TraceField"'+seperator+'"Volume"'+seperator+'"Manifold"'+'\n')
            for p in self.get_polys():
                for v in self.get_volumes(p):
                    for m in self.get_manifolds(p,v):
                        f.write('"'+p+'"'+seperator)
                        f.write('"'+v+'"'+seperator)
                        f.write('"'+m+'"\n')
        finally:
            if type(output_file) == str and f:
                f.close()

    # This filter removes some polynomials with no subfields of degree <= maxsfdegree
    # it doesn't get them all, but does avoid calling nfsubfields; it is quick and approximate.
    def filter_fields(self, maxsfdegree):
        def _filter(p): # for a double break
            deg = pari(p).poldegree()
            for n in xrange(maxsfdegree):
                if _binmiss(n+1,deg):
                    return
            del self.data[p]
        for p in self.data.keys():
            _filter(p)

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
            try:
                data.setdefault(w[0],dict()).setdefault(w[1],list()).append(w[2])
            except IndexError:  # This was initially for debugging, but it can be useful if you grabbed a file while a line was being written 
                print 'Malformed Input: '+str(w)
                continue
        return VolumeData(data = data)
    finally:
        if type(infile) == str and f:
            f.close()
