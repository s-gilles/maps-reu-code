#!/usr/bin/python

from snappy.database import OrientableCuspedCensus, HTLinkExteriors, LinkExteriors
from fractions import gcd

# A simple generator for all manifolds in OrientableCuspedCensus,
# LinkExteriors, HTLinkExteriors
class SimpleGenerator:
    mnum = 0
    pulling_from = None
    pulling_from_str = None

    def __init__(self):
        self.pulling_from = iter(OrientableCuspedCensus)
        self.pulling_from_str = 'ocm'

    def __iter__(self):
        return self

    def next(self):
        while True:
            if self.pulling_from is None:
                raise StopIteration
            try:
                return self.pulling_from.next()
            except:
                if self.pulling_from_str is 'ocm':
                    self.pulling_from = iter(LinkExteriors)
                    self.pulling_from_str = 'le'
                elif self.pulling_from_str is 'le':
                    self.pulling_from = iter(HTLinkExteriors)
                    self.pulling_from_str = 'hle'
                else:
                    self.pulling_from = None
                continue

    def fast_forward_over(self, already_seen_mnum):
        while self.mnum < already_seen_mnum:
            self.mnum = self.mnum + 1
            self.next()

    def next_batch(self, batch_size):
        return [ self.next() ]

# class DehnFillGenerator:
#     all_pqs = None
#     pqs_for_already_seen_range = None
#     already_seen_mnum = 0
#     mnum = 0
#     pulling_from = None
#     pulling_from_str = None
#     current_manifold = None
#     pq_iter = None
# 
#     def pqs_in_range(dehn_pq_limit):
#         pqs = list()
#         for p in range(-1 * DEHN_PQ_LIMIT, DEHN_PQ_LIMIT):
#             for q in range(-1 * DEHN_PQ_LIMIT, DEHN_PQ_LIMIT):
#                 if abs(gcd(p,q)) is 1:
#                     pqs.append((p,q))
#         return pqs
# 
#     def __init__(self, inst_dehn_pq_limit = 24):
#         self.all_pqs = pqs_in_range(inst_dehn_pq_limit)
#         self.pqs_for_already_seen_range = list()
#         self.already_seen_mnum = 0
#         self.mnum = 0
#         self.current_manifold = None
#         self.pulling_from = iter(OrientableCuspedCensus)
#         self.pulling_from_str = 'ocm'
#         pq_iter = 
# 
#     def fast_forward_over(self, inst_already_seen_mnum, inst_dehn_limit):
#         self.already_seen_mnum = inst_already_seen_mnum
#         self.pqs_for_already_seen_range = [ x for x in all_pqs if x not in pqs_in_range(inst_dehn_limit) ]
# 
#     def next(self):
#         while True:
#             try:
#                 
#             if self.pulling_from is None:
#                 raise StopIteration
#             try:
#                 return self.pulling_from.next()
#             except:
#                 if self.pulling_from_str is 'ocm':
#                     self.pulling_from = iter(LinkExteriors)
#                     self.pulling_from_str = 'le'
#                 elif self.pulling_from_str is 'le':
#                     self.pulling_from = iter(HTLinkExteriors)
#                     self.pulling_from_str = 'hle'
#                 else:
#                     self.pulling_from = None
#                 continue
