#!/usr/bin/python

from snappy import *
from snappy.database import OrientableCuspedCensus, HTLinkExteriors, LinkExteriors
from fractions import gcd
from itertools import permutations

# Adds a cache to an iterator that makes it remember the last value sent.
class LastIterator:
    def __init__(self,iterator):
        self.it = iterator
        self.last = None

    def __iter__(self):
        return self

    def last(self):
        return self.last

    def next(self, default = None):
        self.last = self.it.next(default = default)
        return self.last

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

    # Gives all lists of self.length with boolean values with up to self.max_true Trues 
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
                        self.switches = permutations(range(self.length),self.ctrue)
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

class FixedTorusBundleIterator:

    def __init__(self,simplices,start=None):
        self.l = simplices
        self.src = MaskIterator(simplices)
        if start is not None:
# Trouble is, it's not obvious how snappy wants these specified,
# and moreover, some redundancy is certainly present.
# So the only reliable way I can see is to iterate directly.
# Partial implementation left in case this is later resolved.
# Sorry :(
#            try:
#                nm = start.name()
#                # Reverse the string construction that gave this manifold.
#                targ = [False]*(len(nm)-4)
#                if nm[2] is '-':
#                    targ[0] = True
#                elif nm[2] is not '+':
#                    raise ValueError
#                for idx in xrange(len(nm[4:])):
#                    if ch is 'L':
#                        targ[idx+1] =
#                for b in self.src:
            try:
                while self.next() is not start:
                    pass
            except StopIteration:
                print 'Warning: tried to start iterator with non-output '+start.name()+'.'
                self.src.reset()    # Act as if start = None


    def __iter__(self):
        return self

    def next(self, default = None):
        while True:
            try:
                try:
                    bstr = self.src.next()
                except StopIteration:
                    if default is not None:
                        return default
                    else:
                        raise StopIteration
                ostr = 'bo'
                if bstr[0]:
                    ostr += '-R'
                else:
                    ostr += '+R'
                for v in bstr:
                    if v:
                        ostr += 'L'
                    else:
                        ostr += 'R'
                out = Manifold(ostr)
                break
            except (IOError, AttributeError, ValueError): # In the rare case the string was invalid.
                continue
        return out

class TorusBundleIterator:
    def __init__(self, start=2): # start may be int (number of simplices) or manifold
        if isinstance(start, Manifold):
            self.l = len(start.name()) - 3
            self.src = FixedTorusBundleIterator(self.l,start=start)
        elif not isinstance(start, int):
            print 'Warning: tried to start a manifold iterator at '+str(start)+'.'
            start = 2
        if isinstance(start, int):
            self.l = start
            self.src = FixedTorusBundleIterator(start)

    def __iter__(self):
            return self

    def next(self, default = None):
        while True:
            try:
                return self.src.next()
            except StopIteration:
                self.l += 1
                self.src = FixedTorusBundleIterator(self.l)
                continue

class FixedDTIterator:
    # Set start to a manifold this iterator will later output, and the iterator will skip to right after it.
    # That way, you can mark a place to return to more efficently.
    # Or it will be more efficent if I bother to parse mask and abs seperately; won't unless this wastes too much time.
    def __init__(self,crossings,start=None):
        self.n = crossings
        # self.k = components   # Not yet
        m = list()
        for i in xrange(2*crossings):
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
                        self.abs_curr = list(abs_perms.next())
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
            return sub.next()
        except StopIteration:
            self.crossings += 1
            self.sub = FixedDTIterator(crossings)
            return self.sub.next()

# A simple generator for all manifolds in OrientableCuspedCensus,
# LinkExteriors, HTLinkExteriors
class SimpleGenerator:
    mnum = 0
    pulling_from = None
    pulling_from_str = None

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

    def next_batch(self, batch_size):
        ret = list()

        # This is awkward because we want to completely exhaust before raising
        # exception
        ret.append(self.next())
        try:
            for i in range(1, batch_size):
                ret.append(self.next())
        except:
            pass
        return ret

class DehnFillGenerator:
    all_pqs = None
    pqs_for_already_seen_range = None
    already_seen_mnum = 0
    mnum = 0
    pq_iter = None

    def pqs_in_range(dehn_pq_limit):
        pqs = list()
        for p in range(-1 * dehn_pq_limit, dehn_pq_limit):
            for q in range(0, dehn_pq_limit):
                if abs(gcd(p,q)) is 1:
                    pqs.append((p,q))
        return pqs

    def __init__(self, full_dehn_pq_limit = 24, already_done_manifolds_up_to = 0):
        self.all_pqs = pqs_in_range(full_dehn_pq_limit)
        pqlen = len(self.all_pqs)
        self.mnum = 0
        self.pulling_from = SimpleGenerator()
        for i in range(0, already_done_manifolds_up_to):
            self.pulling_from.next()
        self.current_manifold = self.pulling_from.next()
        pq_iter = iter(all_pqs)

    def next(self):
        while True:
            try:
                p, q = pq_iter.next()
                m = self.current_manifold.copy()
                m.dehn_fill((p,q))
                return m
            except StopIteration:
                try:
                    self.current_manifold = self.pulling_from.next()
                    pq_iter = iter(all_pqs)
                except StopIteration:
                    raise StopIteration
