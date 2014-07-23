#!/usr/bin/python

from snappy import *
from snappy.database import OrientableCuspedCensus, HTLinkExteriors, LinkExteriors
from fractions import gcd
from itertools import permutations, combinations, product

# Default maximum number of strands in a braid; high numbers of strands tend to make snappy take forever.
DEF_STRANDS = 4

def get_block(blockidx, blocksz, obj):
    """
    Given obj, a list, return the intersection of
    obj[blockidx*blocksz:(blockidx+1)*blocksz] and obj

    Ex: get_block(2, 100, range(250) returns [200, 201, ..., 249]
    """
    if blockidx*blocksz > len(obj):
        return []
    elif (blockidx+1)*blocksz > len(obj):
        return obj[blockidx*blocksz:]
    else:
        return obj[blockidx*blocksz:(blockidx+1)*blocksz]

class ForwardingIterator:
    """
    An iterator that, given a source iterator and a function, returns
    the result of applying the function to each element of the source
    iterator in turn.

    Ex: ForwardingIterator(iter(range(10)), lambda x:x**2) will return
    0, 1, 4, 9, ... 81
    """
    def __init__(self,source,funct):
        self.source = source
        self.funct = funct

    def __iter__(self):
        return self

    def next(self):
        return self.funct(self.source.next())

class StartIterator:
    """
    An iterator that, given a source iterator and a starting point,
    skips forward until the result of next() would be the input
    provided.

    Note: the output and next() elements are compared using the ==
    operator.

    Ex: StartIterator(iter(range(10)), 3) will return 3, 4, ... 9
    """
    def __init__(self,source,output):
        self.source = source
        self.return_special = True
        self.special_first_element = output
        try:
            while True:
                if source.next() == output:
                    break
        except StopIteration:
            pass    # This way we only throw StopIteration when next is called.

    def __iter__(self):
        return self

    def next(self):
        if self.return_special:
            self.return_special = False
            return self.special_first_element
        return self.source.next()

class BatchIterator:
    """A wrapper around an iterator (source) that adds a next_batch()
    method. This method clumps together results into lists of size
    batch_size and returns them.

    Note that if a batch would return more values than are left in the
    iterator, but at least one value would be returned, the remainder of
    the iterator is returned instead of raising a StopIteration

    Ex: BatchIterator(iter(range(50)), 15) will return (using
    next_batch() instead of next()):

    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
    [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
    [30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44],
    [45, 46, 47, 48, 49]
    """
    def __init__(self,source,batch_size):
        self.source = source
        self.def_batch_size = batch_size

    def next(self):
            return self.source.next()

    def next_batch(self, batch_size = -1):  # Default value that indicates preset batch_size.
        """
        Return the next batch.  Optionally, the parameter batch_size can
        be provided, which, if positive, overrides the batch_size that
        was provided in the initializer for this method.

        Note: If the source iterator cannot provide the requested number
        of elements, the entire remainder of the iterator will be
        returned. Only if the iterator has no elements left will a
        StopIteration be raised.
        """
        if batch_size <= 0:
            batch_size = self.def_batch_size
        ret = [self.source.next()]   # Raise StopIteration if nothing to return at all.
        for i in xrange(batch_size - 1):
            try:
                ret.append(self.source.next())
            except StopIteration:
                break   #...but if we have a partial batch, return it.
        return ret

class MaskIterator:
    """
    Outputs all possible bitmasks of a given length, given as arrays of
    boolean values

    Ex: MaskIterator(3) returns [False, False, False], [True, False,
    False], [False, True, False], ...,  [True, True, True]
    """
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
    def next(self):
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
                        raise StopIteration
            return out

    def reset(self):
        """
        Reset this iterator to the beginning of its run
        """
        self.ctrue = 0
        self.switches = iter([[]])

def _get_bool_array(integer, length):
    """
    Returns true/false list representation of the integer with desired
    length (padding falses)
    """
    s = bin(integer)[2:]
    if length < len(s): # Integer is too big
        raise ValueError
    o = [False]*length
    for i in xrange(len(s)):
        if s[-(i+1)] == '1':
            o[-(i+1)] = True
    return o

def _something_new(bool_array):
    """
    Returns False if the array is a duplicate of a smaller sub-array,
    True if it cannot be constructed as a duplicate of any of its
    sub-components.

    Ex: _something_new([False, False, False, True]) is True, but
    _something_new([False, True, False, True, False, Treu]) is False,
    since the latter is [Fale, True] * 3.
    """
    x = list(bool_array)
    for k in xrange(1, len(x)):
        if len(x) % k == 0:
            if x == x[:k] * (len(x)/k):
                return False
    return True

class RandomIterator:
    """
    Accepts a source, which must be an iterator that returns manifolds,
    and an optional number max_tries (defaulting to 16). This
    transparently returns manifolds from the source iterator, but each
    manifold that would be returned as not geometric is first randomized
    (by the Manifold.randomize() method) for up to max_tries attempts or
    until the manifold becomes geometric.
    """
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
        """
        Return a list of all manifolds that this iterator has returned
        that could not be randomized into a geometric manifold.
        """
        return self.failures


class FixedTorusBundleIterator:
    """
    Given the following:

    start_index: a positive integer equal to the input for
    _get_bool_array to get the desired manifold. (valid values are in
    [0, 2**simplices) )

    skip: an optional boolean, defaulting to False.  If set to True, the
    iterator will skip manifolds that are multiples of others (as
    defined by a false result to _something_new)
    """
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
                    bstr = _get_bool_array(self.idx, self.l)
                    self.idx += 1
                    if self.skip:   # try some labor saving
                        while not _something_new([False]+bstr[1:]):
                            bstr = _get_bool_array(self.idx, self.l)
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

    def last_idx(self):
        """
        Return the last index that this iterator computed
        """
        return self.idx - 1  # idx gets incremented in next() call before it is used

class TorusBundleIterator:
    """
    Given a start_length (defaulting to 2), start_idx (defaulting to 0),
    and skip (defaulting to False), return Torus Bundle-based manifolds.
    The Torus Bundles start as defined with codes of length
    start_length, and within that start_length they start at the
    position given by start_idx, and if skip is true, Torus Bundles that
    are known to be multiples of other Torus Bundles are skipped.
    """
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
    # Now, invert _get_bool_array
    for n in xrange(len(bits)):
        if bits[-(n+1)]:
            idx += 2**n
    return idx
    "b+-RLRRRRLLRRL(1,1)"

class FixedDTIterator:
    """
    Returns all DT iterators that have a fixed number of
    crossings. Given an optional start element, skip forward until that
    element would be returned.

    Note: This iterator has issues with SnapPy's usage of SQLite3.  It
    should be considered alpha quality, and not actually used.
    """
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
    """
    Return all DT codes, starting with three crossings.

    Note: This iterator has issues with SnapPy's usage of SQLite3.  It
    should be considered alpha quality, and not actually used.
    """
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

class FixedBraidIterator:
    """
    Generates manifolds from braid words of a given word length and a
    number of strands.

    Note: These are constructed as passing to SnapPy the string 'Braid:'
    followed by a string of integers where x is never followed by -x,
    and the first digit, in this iterator, is restricted to always be
    1. This is conjectured to match all Braid-generated manifolds by
    isomorphism.
    """
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
    """
    Starting at length start_length (defaulting to to), return all braids.

    Note: These are constructed as passing to SnapPy the string
    'Braid:' followed by a string of integers where x is never followed
    by -x, and the first digit, in this iterator, is restricted to
    always be 1. This is conjectured to match all Braid-generated
    manifolds by isomorphism.
    """
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

def get_braid_idx(man_nm, strands = DEF_STRANDS):
    """
    Return the index used to generate a given braid with the given
    number of strands.
    """
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

def _pqs_in_range(dehn_pq_limit, num_cusps):
    """
    Return an iterator.  This iterator, at each step, returns a
    tuple. The contents of this tuple are num_cusps other tuples, and
    each of these is of the form (p,q), where 0 <= p <= dehn_pq_limit,
    -dehn_pq_limit <= q <= dehn_pq_limit, and gcd(p,q) <= 1.

    Ex: pqs_in_range(3, 2) returns
    ((-3, 1), (-3, 1)),
    ((-3, 1), (-3, 2)),
    ((-3, 1), (-2, 1)),
    ((-3, 1), (-2, 3)),
    ((-3, 1), (-1, 0)),
    ...
    ((3, 2), (2, 1))
    ((3, 2), (2, 3))
    ((3, 2), (3, 1))
    ((3, 2), (3, 2))
    """
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
    """
    Given a source that returns manifold objects, return all reasonable
    dehn fillings of those manifolds.

    full_dehn_pq_limit is a list, such that full_dehn_pq_limit[i] is the
    pq_limit of a manifold with i cusps (defaulting to dehn_pq_limit[0]
    if i would cause an IndexException)

    For example, DehnFillIterator(iter([Manifold('m004'),
    Manifold('s776')])) will return:

    m005(-16,1),
    m004(-16,3)
    m004(-16,5)
    m004(-16,7)
    m004(-16,9)
    ...
    m004(16,13)
    m004(16,15)
    s776(-8,1)(-8,1)(-8,1)
    s776(-8,1)(-8,1)(-8,3)
    s776(-8,1)(-8,1)(-8,5)
    ...
    s776(8,1)(5,7)(-2,1)
    s776(8,1)(5,7)(-2,3)
    s776(8,1)(5,7)(-2,5)
    s776(8,1)(5,7)(-2,7)
    s776(8,1)(5,7)(-1,0)
    ...

    This will exhaust all dehn fillings of all manifolds, with each dehn
    filling (p,q) where p is constrained by [0, full_dehn_pq_limit[i]]
    and q is constrained by [-full_dehn_pq_limit[i],
    full_dehn_pq_limit[i]], and gcd(p,q) <= 1.

    Note that (0,0) is included as a dehn filling, which will return a
    manifold equal to the unsurgeried dehn filling.
    """
    def __init__(self, source, full_dehn_pq_limit = [6, 16, 12, 8, 6, 4, 3, 3, 2, 2, 2], fast_forward_to_pq = None):
        self.mnum = 0
        self.pulling_from = source
        self.current_manifold = self.pulling_from.next()
        try:
            pq_limit = full_dehn_pq_limit[self.current_manifold.num_cusps()]
        except:
            pq_limit = full_dehn_pq_limit[0]
        self.all_pqs = _pqs_in_range(pq_limit, self.current_manifold.num_cusps())
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
                self.all_pqs = _pqs_in_range(pq_limit, self.current_manifold.num_cusps())
                self.pq_iter = iter(self.all_pqs)
