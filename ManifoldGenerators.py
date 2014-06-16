#!/usr/bin/python

from snappy.database import OrientableCuspedCensus, HTLinkExteriors, LinkExteriors
from fractions import gcd

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
