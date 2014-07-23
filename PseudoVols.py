#!/usr/bin/python

from snappy import Manifold
from cypari import *
from ManifoldIterators import *
from multiprocessing import *
import copy

# Until this gets globalized
EPSILON = .0000000000001    # same as VolumeProcessing
MAX_ITF = 8                 # ^
pari.set_real_precision(100)

def prepare_pvolume_file(maniter, ofilenm, append = False, engine = 'magma', max_secs = 20, retrieve = True, period = 100, seperator = ';'):
    """The same as calling get_volume_data(mans).write_to_csv(ofilenm) with the given parameters,
except output will be written out every period manifolds and logs generated, instead of all at once."""

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
                v = get_volume_data(ForwardingIterator(block.__iter__(),lambda m : str(m)),engine=engine,max_secs=max_secs,retrieve=retrieve)
                v.write_to_csv(f,append=append,seperator=seperator)
                append = True  # we must be appending after the first time
                ctr = 0
                block = list()
                if done:
                    break
    finally:
        f.close()

# Count the number of distinct non-zero (<EPSILON) values up to sign
def _distinct_abs(vol_list):
    return len(set([abs(v) for v in vol_list if abs(v) >= EPSILON]))   # remove zero volumes

def get_volume_data(man_nms, engine = 'magma', max_secs = 20, retrieve = True, max_itf_degree = MAX_ITF):

    """ Returns a VolumeData object containing exotic volumes for manifolds with the given names
Volumes' precision is based on pari, so set it there
set retrieve = False to skip retrieving ptolemy data from files
set engine = None to skip computing ptolemy data in an engine
set max_secs to specify how long we will spend computing a given manifolds' data before killing the engine and moving on;
  specifying None means we will never give up (unless something crashes)
if the engine given crashes, so will IDLE and SnapPy; to avoid this, run this command only from within python scripts.
Manifolds with more than floor(max_itf_degree/2) distinct volumes to an obstruction class
will have their data for that obstruction class removed, since this demonstrates an invariant trace field with too high ncp
Set to None and it will be ignored.""" 

#TODO: special case max_secs=None to not bother with processes

    if engine:
        def _use_engine(v,p):   # this function will be called in a second process to facilitate time limits
            p.send(v.compute_solutions(engine = engine))
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
                    mine, theirs = Pipe(duplex = False)
                    p = Process(target=_use_engine,args=[var,theirs])
                    p.daemon = True
                    p.start()
                    if mine.poll(max_secs): # Here is the time limit stuff
                        sols = mine.recv()
                        p.terminate()
                    else:
                        p.terminate()   # give up on this one
                        print 'Computation took too long; skipping '+nm
                        continue
                else:
                    print 'No engine and no data retrieved; skipping '+nm
                    continue
            if sols:
                pairs = [(c.number_field(),c.volume_numerical()) for c in sols.flatten()]   # TODO does this flatten do anything?
                data = list()
                for p in pairs:
                    if _distinct_abs(p[1]) <= max_itf_degree/2:
                        data.append(p)
                    else:
                        data.append(None)
                    # data = [p for p in data if _distinct_abs(p[1]) <= max_itf_degree/2] # arc
                for cl_idx in xrange(len(data)):
                    if data[cl_idx]:
                        for v in data[cl_idx][1]:
                            recs.setdefault(str(data[cl_idx][0]),dict()).setdefault(str(v),list()).append((nm,cl_idx))
            else:
                print 'Got no solutions; skipping '+nm
        except Exception as e:
            print(str(e))+'; skipping '+nm
            continue
    for p in recs.keys():
        for v in recs[p].keys():
            recs[p][v] = list(set(recs[p][v]))
    return VolumeData(data = recs)

def get_potential_trace_fields(poly):
    """Given a minimal polynomial of a trace field, returns a list of minimal polynomials of the potential invariant trace fields."""
    pol = pari(poly)
    try:
        return [str(rec[0].polredabs()) for rec in pol.nfsubfields()[1:] if _binmiss(rec[0].poldegree(),pol.poldegree())]   # poldegree returns int
    except: # we want cypari.gen.PariError, but no idea how to reference; fortunately, anything else will just raise again
        pol = pol.polredabs()
        return get_potential_trace_fields(str(pol))

def _binmiss(s,l):
    while s < l:
        s *= 2
    return s == l

# Wrapper for manipulating data on pseudo-volumes
class VolumeData:
    """This class is for storage and manipulation of exotic volumes of some manifolds.
Given a value for data, the constructor makes a VolumeData object wrapping it.
The datastructure is {poly:{volume:(manifold,obstruction_class_index)}}
It's usually not nescecary to make these yourself; collection and read methods return them for you."""

    # structure: dict poly ---> dict volume ---> [(manifold,obstr_cl)]
    def __init__(self, data = dict()):
        self.data = data

    def get_polys(self):
        """Returns (as a list of strings) the minimal polynomials of the ptolemy/trace fields for the volumes in this object."""
        return self.data.keys()

    def get_volumes(self,poly):
        """Returns (as a list of strings) the volumes that occur over the field with the given minimal polynomial."""
        return self.data[poly].keys()

    def get_manifolds(self,poly,volume):
        """Returns a list of the names of manifolds that produce the given minimal polynomial/volume pair.""" 
        return [p[0] for p in self.data[poly][volume]]

    def combine_with(self,other):
        """Returns a VolumeData object containing the data from self and other; in case of a conflict (which should not occur), 
        the other's data takes precedence."""
        new_data = copy.deepcopy(self.data)
        for p in other.get_polys():
            for v in other.get_volumes(p):
                new_data.setdefault(p,dict()).setdefault(v,list()).extend(other.data[p][v])
        return VolumeData(data = new_data)

    # given an (open, ready to write) file object or valid filename, writes the data
    def write_to_csv(self, output_file, seperator = ';', append = False):
        """Writes out the data to output_file, provided output_file is a valid (open, ready to write) File object or filename."""
        try:
            if type(output_file) == str:
                if append:
                    f = open(output_file,'a')
                else:
                    f = open(output_file,'w')
            else:
                f = output_file
            if not append:
                f.write('"TraceField"'+seperator+'"Volume"'+seperator+'"Manifold"'+seperator+'"ObstructionClass"'+'\n')
            for p in self.get_polys():
                for v in self.get_volumes(p):
                    for rec in self.data[p][v]:
                        f.write('"'+p+'"'+seperator)
                        f.write('"'+v+'"'+seperator)
                        f.write('"'+rec[0]+'"'+seperator)
                        f.write('"'+str(rec[1])+'"\n')
        finally:
            if type(output_file) == str and f:
                f.close()

    def filter_fields(self, maxsfdegree=MAX_ITF):
        """This filter removes some polynomials with no subfields of degree <= maxsfdegree
        it doesn't get them all, but does avoid calling nfsubfields; it is quick and approximate."""
        def _filter(p): # for a double break
            deg = pari(p).poldegree()
            for n in xrange(maxsfdegree):
                if _binmiss(n+1,deg):
                    return
            del self.data[p]
        for p in self.data.keys():
            _filter(p)
    
    # Remove all volumes that are integral multiples of another volume (including 1*)
    # To register as an integral multiple, the decimal part of big/small must be less than epsilon
    # Will remove manifolds if all their pvols were integral multiples of other pvols
    def cull_volumes(self, epsilon = EPSILON):  # code adapted from VolumeProcessing
        for poly in self.get_polys():        
            vols = self.get_volumes(poly)
            i = 0
            while i < len(vols) - 1:
                j = i + 1
                while j < len(vols):
                    try:
                        if is_int(float(vols[i])/float(vols[j]), epsilon = epsilon) and gen.pari(vols[i] + ' > ' + vols[j]) == 1:
                            # We have to throw away (culled) manifold names to let all culled manifolds have the same volume
                            # [j] divides [i] so remove [i]
                            del self.data[poly][vols.pop(i)]
                            # i is already effectivley incremented, so we must offset it
                            i = i-1
                            break
                        elif is_int(float(vols[j])/float(vols[i]), epsilon = epsilon):
                            # this time, remove [j]
                            del self.data[poly][vols.pop(i)]
                            # j is effectivley incremented, no need to do it
                        else:
                            j += 1
                    except (ValueError, ZeroDivisionError): # bad quotient; not a linear combination either way so...
                        j += 1
                i += 1

    def remove_nonpositive_vols(self, epsilon = EPSILON):
        """Removes any volume less than epsilon"""
        for p in self.get_polys():
            for v in self.get_volumes(p):
                try:
                    if float(v) < epsilon:
                        print str(v)+'was less than '+str(epsilon)
                        del self.data[p][v]
                except: # v was really close to 0
                    del self.data[p][v]
                    print str(v)+'threw a value error: '+str(e)

    
    def clean(self, maxsfdegree=MAX_ITF, epsilon = EPSILON):
        """Runs several methods for decreasing size without losing much information
        Set maxsfdegree to None to avoid culling based on subfield degree."""
        if maxsfdegree:
            self.filter_fields(maxsfdegree)
        self.remove_nonpositive_vols(epsilon = epsilon)

    # Cut down to 1 manifold per poly,vol pair.
    def remove_duplicate_manifolds(self):
        for p in self.get_polys():
            for v in self.get_volumes(p):
                self.data[p][v] = [self.data[p][v][0]]

def is_int(fl, epsilon = EPSILON):
    return fl % 1 < epsilon or 1 - (fl % 1) < epsilon

def read_volumedata_csv(infile, seperator = ';'):
    """Given an (open, ready to read data) file object or valid filename, reads the file and returns a VolumeData that would write it."""
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
                data.setdefault(w[0],dict()).setdefault(w[1],list()).append((w[2],int(w[3])))
            except IndexError:  # This was initially for debugging, but it can be useful if you grabbed a file while a line was being written 
                print 'Malformed Input: '+str(w)
                continue
        return VolumeData(data = data)
    finally:
        if type(infile) == str and f:
            f.close()
