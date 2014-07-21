#!/usr/bin/python

from snappy import *
from snappy.database import OrientableCuspedCensus, HTLinkExteriors, LinkExteriors
from fractions import gcd
from itertools import permutations, combinations, product

# Default maximum number of strands in a braid; high numbers of strands tend to make snappy take forever.
DEF_STRANDS = 4

# Given a list-like obj, returns the intersection of obj[blockidx*blocksz:(blockidx+1)*blocksz] and obj;
# i.e. if blocksz is 100, you get obj[0:100],obj[100:200],obj[200,300]... for increasing blockidx
def get_block(blockidx,blocksz,obj):
    if blockidx*blocksz > len(obj):
        return []
    elif (blockidx+1)*blocksz > len(obj):
        return obj[blockidx*blocksz:]
    else:
        return obj[blockidx*blocksz:(blockidx+1)*blocksz]

# This iterator applies a given funct to source.next() before returning it.
class ForwardingIterator:
    def __init__(self,source,funct):
        self.source = source

    def __iter__(self):
        return self

    def next(self):
        return str(self.source.next())

# Skips the source iterator to after the given output
class StartIterator:
    def __init__(self,source,output):
        self.source = source
        try:
            while True:
                if source.next() == output:
                    break
        except StopIteration:
            pass    # This way we only throw StopIteration when next is called.

    def __iter__(self):
        return self

    def next(self):
        return self.source.next()

# Iterates through each element of 'sources' in order.
class QueueIterator:
    def __init__(self,sources):
        self.sources = sources
        self.finished = list()

    def __iter__(self):
        return self

    def next(self):
        while True:
            try:
                return self.sources[0].next()
            except StopIteration:
                self.finished.append(self.sources.pop(0))
                continue
            except IndexError: # Out of stuff to return
                raise StopIteration

    # Adds the source at the given index, bumping everything back
    def add(source, index):
        if index < 0:
            index = 0
        elif index > len(self.sources):
            index = len(self.sources)
        self.sources.insert(source,index)

    # Removes and returns a source
    def pop(index):
        try:
            return self.sources.pop(index)
        except IndexError:
            pass

    # You could modify the list, but please don't.
    def peek():
        return self.sources

    # Lets you see the current iterator
    def current():
        return self.sources[0]

    # A list of iterators we went through.
    def finished():
        return self.finished

# Adds a cache to an iterator that makes it remember the last value sent.
class LastIterator:
    def __init__(self,iterator):
        self.it = iterator
        self.last_item = None

    def __iter__(self):
        return self

    def last(self):
        return self.last_item

    def next(self, default = None):
        try:
            self.last_item = self.it.next()
        except StopIteration:
            if default is not None:
                return default
            else:
                raise
        return self.last_item

# Generalized ability to grab a batch from an iterator 'source'.
class BatchIterator:
    def __init__(self,source,batch_size):
        self.source = source
        self.def_batch_size = batch_size

    def next(self):
            return self.source.next()

    def next_batch(self, batch_size = -1):  # Default value that indicates preset batch_size.
        if batch_size <= 0:
            batch_size = self.def_batch_size
        ret = [self.source.next()]   # Raise StopIteration if nothing to return at all.
        for i in xrange(batch_size - 1):
            try:
                ret.append(self.source.next())
            except StopIteration:
                break   #...but if we have a partial batch, return it.
        return ret

    def last(self):
        return self.source.last()

# Because itertools returns the same thing more than once.
class MaskIterator:

    def __init__(self,length, max_true = None):
        self.length = length
        self.max_true = max_true
        if max_true == None:
            self.max_true = length
        self.ctrue = 0
        self.switches = iter([[]])

    def __iter__(self):
        return self

    # Gives all lists of self.length with boolean values with up to
    # self.max_true Trues
    def next(self, default = None):
            # produce output
            out = [False]*self.length
            while True:
                try:
                    for i in self.switches.next():
                        out[i] = True
                    break
                except StopIteration:
                    self.ctrue += 1
                    if self.ctrue <= self.max_true:
                        self.switches = combinations(range(self.length),self.ctrue)
                        continue
                    else:
                        if default is None:
                            raise StopIteration
                        else:
                            return default
            return out

    def reset(self):
        self.ctrue = 0
        self.switches = iter([[]])

# Returns true/false list representation of the integer with desired length (padding falses)
def get_bool_array(integer, length):
    s = bin(integer)[2:]
    if length < len(s): # Integer is too big
        raise ValueError
    o = [False]*length
    for i in xrange(len(s)):
        if s[-(i+1)] == '1':
            o[-(i+1)] = True
    return o

def something_new(bool_array):
    x = list(bool_array)
    for k in xrange(1, len(x)):
        if len(x) % k == 0:
            if x == x[:k] * (len(x)/k):
                return False
    return True

# Randomizes each manifold in man_iter up to max_tries times to get a geometric solution, then gives up and moves on 
class RandomIterator:
    def __init__(self, man_iter, max_tries = 16):  
        self.source = man_iter
        self.max_tries = max_tries
        self.failures = list()

    def __iter__(self):
        return self

    def next(self):
        man = self.source.next()
        oman = man.copy()
        tries = 0
        while oman.solution_type(enum = True) != 1:
            oman.randomize()
            tries += 1
            if tries == self.max_tries:
                self.failures.append(man)
                man = self.source.next()
                oman = man.copy()
                tries = 0
        return oman

    # Gives a list of the manifolds that we couldn't make geometric.
    def get_failures(self):
        return self.failures
        

class FixedTorusBundleIterator:
    # start_index is a positive integer equal to the input for get_bool_array to get the desired manifold
    # valid values are 0 to 2**simplices - 1
    # (probably actually 1 to 2**simplices - 2, but the extreme values will just be skipped naturally) 
    # If skip = True, the iterator will skip integers that produce a multiple of some other volume with same or lesser word length.
    # Testing for this may be expensive (who knows?), but cost probably grows slower than volume computation's.
    def __init__(self, simplices, start_index=0, skip=False):
        self.l = simplices
        self.idx = start_index
        self.skip = skip

    def __iter__(self):
        return self

    def next(self, default = None):
        while True:
            try:
                try:
                    bstr = get_bool_array(self.idx, self.l)
                    self.idx += 1
                    if self.skip:   # try some labor saving
                        while not something_new([False]+bstr[1:]):
                            bstr = get_bool_array(self.idx, self.l)
                            self.idx += 1
                except ValueError:  # We must be done since idx >= 2**simplices
                    if default is not None:
                        return default
                    else:
                        raise StopIteration
                ostr = 'bo'
                if bstr[0]:
                    ostr += '-R'
                else:
                    ostr += '+R'
                for v in bstr[1:]:
                    if v:
                        ostr += 'L'
                    else:
                        ostr += 'R'
                out = Manifold(ostr)
                break
            except (IOError, AttributeError, ValueError): # In the rare case the string was invalid.
                continue
        return out


    # Note that this method does not account for the case where the last index was skipped because it produced a ValueError
    # This seems to happen whenever that index = 2**n or 2**n + 1
    # Note also that this and other last_foo methods are unreliable when next() has never been called.
    def last_idx(self):
        return self.idx - 1  # idx gets incremented in next() call before it is used

class TorusBundleIterator:
    def __init__(self, start_length=2, start_idx=0, skip = False):
        self.l = start_length
        self.skip = skip
        self.src = FixedTorusBundleIterator(self.l, start_index = start_idx, skip = self.skip)

    def __iter__(self):
            return self

    def next(self, default = None):
        while True:
            try:
                return self.src.next()
            except StopIteration:
                self.l += 1
                self.src = FixedTorusBundleIterator(self.l, skip = self.skip)
                continue

    def last_idx(self):
        i = self.src.last_idx()
        if i == -1:
            return 2**(self.l-1)-1 # We just got to this length
        else:
            return i

    def last_length(self):
        if self.src.last_idx() == -1:
            return self.l - 1
        else:
            return self.l

def get_torus_idx(man_nm):
    nm = man_nm[2:] # assume that nothing before the sign matters
    # Invert next (in the fixed iterator)
    bits = list()
    if nm[0] == '-':
        bits += [True]
    else:
        bits += [False]
    for x in nm[2:]:    # skip the fixed R, it makes no bit
        if x == 'L':
            bits += [True]
        else:
            bits += [False]
    idx = 0
    # Now, invert get_bool_array
    for n in xrange(len(bits)):
        if bits[-(n+1)]:
            idx += 2**n
    return idx
    "b+-RLRRRRLLRRL(1,1)"

class FixedDTIterator:
    # Set start to a manifold this iterator will later output, and the iterator
    # will skip to right after it.  That way, you can mark a place to return to
    # more efficently.  Or it will be more efficent if I bother to parse mask
    # and abs seperately; won't unless this wastes too much time.
    def __init__(self,crossings,start=None):
        self.n = crossings
        # self.k = components   # Not yet
        m = list()
        for i in xrange(crossings):
            m.append(2*(i))
        self.abs_perms = permutations(m,self.n)
        self.abs_curr = list(self.abs_perms.next())
        self.sgn_perms = MaskIterator(self.n,self.n//2)  # mirrors don't count
        # In case they want us to start:
        if start is not None:
            try:
                while start is not self.next():
                    pass
            except StopIteration:
                self.abs_perms = permutations(m,self.n)
                self.sgn_perms.reset()
        self.abs_curr = list(self.abs_perms.next())

    def __iter__(self):
        return self

    def next(self, default = None):
        while True:
            try:
                try:
                    self.sgn_curr = list(self.sgn_perms.next())
                except StopIteration:
                    try:
                        self.abs_curr = list(self.abs_perms.next())
                    except StopIteration:
                        if default is None:
                            raise StopIteration
                        else:
                            return default
                    self.sgn_perms.reset()
                    self.sgn_curr = self.sgn_perms.next()
                for i in xrange(len(self.abs_curr)):
                    if self.sgn_curr[i]:
                        self.abs_curr[i] = self.abs_curr[i] * -1
                print 'DT:['+str(self.abs_curr)+']' # DEBUG
                return Manifold('DT:[('+str(self.abs_curr).replace(' ','').strip('[]')+')]')  # HERE fix formatting
                # Looks like this fails randomly according to unexplained constraints
                # Chuck in a loop that catches relevant exceptions...
            except AttributeError, ValueError:
                print 'Failed: \'DT:[('+str(self.abs_curr).replace(' ','').strip('[]')+')]\'.'
                continue

class DTIterator:
    def __init__(self,start=None):
        self.crossings = 3
        self.sub = FixedDTIterator(self.crossings)
        # Implement start mattering

    def __iter__(self):
        return self

    def next(self, default = None):
        try:
            return self.sub.next()
        except StopIteration:
            self.crossings += 1
            self.sub = FixedDTIterator(crossings)
            return self.sub.next()

# Generates manifolds from braid words of a given word length and number of strands
class FixedBraidIterator:

    # Large values for strands can make computations take a very long time
    def __init__(self, word_length, strands = DEF_STRANDS, start_idx = 0):
        self.length = word_length
        self.base = 2*strands - 1
        self.max_pow = word_length - 2 
        self.idx = start_idx
        self.strands = strands
        self.stop_idx = self.base**(self.max_pow + 1) # only need to compute this once

    def __iter__(self):
        return self

    def next(self):
        if self.idx >= self.stop_idx:
            raise StopIteration
        word = [1]
        power = self.max_pow
        rem = self.idx
        while power >= 0:
            place = self.base**power
            word.append(rem//place) # get the next digit
            rem = rem % place       # get ready to recurse
            power -= 1
        for n in xrange(1,len(word)):
            word[n] -= self.strands     # make range -n to n-2 instead of 0 to 2n-1
            if word[n] >= 0:            # now -n to -1 and 1 to n-1 
                word[n] += 1
            if word[n] >= -1*word[n-1]: # bump up to avoid ...n,-n... which would cancel out; now -n to -1 and 1 to n
                word[n] += 1
            if word[n] == 0:            # what if word[n-1] >= 1 and we just sent -1 to 0?
                word[n] = 1
        self.idx += 1    # increment after we work
        return Manifold('Braid:'+str(word).replace(' ',''))
        
# Generates manifolds from braids with a given number of strands
class BraidIterator:

    # Large values for strands can make computations take a very long time
    def __init__(self, start_length=2, strands = DEF_STRANDS, start_idx=0):
        self.length = start_length
        self.source = FixedBraidIterator(self.length, start_idx = start_idx)
        self.strands = strands

    def __iter__(self):
        return self

    def next(self):
        try:
            return self.source.next()
        except StopIteration:
            self.length += 1
            self.source = FixedBraidIterator(self.length, strands = self.strands)
            return self.next()

# Returns the index used to generate a given braid with the given number of strands.
def get_braid_idx(man_nm, strands = DEF_STRANDS):
    if man_nm[:6] != 'Braid:':  # make sure input is valid
        raise ValueError
    digits = [int(x) for x in man_nm[6:].replace(' ','').strip('[]').split(',')]    # cut off first digit
    orig = list()
    for n in xrange(1,len(digits)):
        res = digits[n]
        if res > -1*digits[n-1]:  # reverse bumping up to avoid -last
            res -= 1
            if res == 0:
                res -= 1
        if res >= 1:  # reverse bumping up to avoid 0
            res -= 1
        res += strands  # reverse shifting to -strands , strands-2
        orig.append(res)
    base = 2*strands - 1    
    idx = 0
    for n in xrange(len(orig)):
        idx += orig[-(n+1)] * (base**n)
    return idx
    

# A simple generator for all manifolds in OrientableCuspedCensus,
# LinkExteriors, HTLinkExteriors
class SimpleIterator:
    def __init__(self, already_seen_mnum = 0):
        self.mnum = 0
        self.pulling_from = iter(OrientableCuspedCensus)
        self.pulling_from_str = 'ocm'
        while self.mnum < already_seen_mnum:
            next()

    def __iter__(self):
        return self

    def next(self):
        while True:
            if self.pulling_from is None:
                raise StopIteration
            try:
                return self.pulling_from.next()
            except StopIteration:
                if self.pulling_from_str is 'ocm':
                    self.pulling_from = iter(LinkExteriors)
                    self.pulling_from_str = 'le'
                elif self.pulling_from_str is 'le':
                    self.pulling_from = iter(HTLinkExteriors)
                    self.pulling_from_str = 'hle'
                else:
                    self.pulling_from = None
                continue

def pqs_in_range(dehn_pq_limit, num_cusps):

    # pqs = [ None, (0,0), (0,1), (0,3), (5,12), (9,13), ... ]
    pqs = list()
    for p in range(-1 * dehn_pq_limit, dehn_pq_limit + 1):
        for q in range(0, dehn_pq_limit + 1):
            if abs(gcd(p,q)) <= 1:
                pqs.append((p,q))

    # pqs_mult = [ pqs, pqs, pqs... ]
    # because we wish to return pqs x pqs x pqs ... x pqs
    pqs_mult = list()
    for i in range(0, num_cusps):
        pqs_mult.append(pqs)
    return product(*pqs_mult)

class DehnFillIterator:
    def __init__(self, source, full_dehn_pq_limit = [6, 16, 12, 8, 6, 4, 3, 3, 2, 2, 2], fast_forward_to_pq = None):
        self.mnum = 0
        self.pulling_from = source
        self.current_manifold = self.pulling_from.next()
        try:
            pq_limit = full_dehn_pq_limit[self.current_manifold.num_cusps()]
        except:
            pq_limit = full_dehn_pq_limit[0]
        self.all_pqs = pqs_in_range(pq_limit, self.current_manifold.num_cusps())
        self.pq_iter = iter(self.all_pqs)
        if fast_forward_to_pq is not None:
            tmp_pq = None
            while tmp_pq != fast_forward_to_pq:
                tmp_pq = self.pq_iter.next()
        self.dehn_pq_limit = full_dehn_pq_limit
        self.peek_ret = None

    def next(self):
        while True:
            try:
                pqs = self.pq_iter.next()
                m = self.current_manifold.copy()

                curr_idx = 0
                for curr_pq in pqs:
                    if curr_pq is not None:
                        m.dehn_fill(curr_pq, curr_idx)
                    curr_idx = curr_idx + 1
                self.peek_ret = m
                return m
            except StopIteration:
                self.current_manifold = self.pulling_from.next()
                try:
                    pq_limit = self.dehn_pq_limit[self.current_manifold.num_cusps()]
                except:
                    pq_limit = self.dehn_pq_limit[0]
                self.all_pqs = pqs_in_range(pq_limit, self.current_manifold.num_cusps())
                self.pq_iter = iter(self.all_pqs)

    def peek(self):
        return self.peek_ret
