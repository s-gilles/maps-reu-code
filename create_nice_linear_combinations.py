#!/usr/bin/python2
import os
from VolumeProcessing import *
from PseudoVols import *

n=2
maxdegree=18
PSEUDO_VOLUME_FILENAME='pseudo_volumes_n='+str(n)+'_nice.csv'

def all_fields(polynomial,trace_field,n):
    return True

# Create pseudo-volume file
if not os.path.exists(PSEUDO_VOLUME_FILENAME):
    prepare_pvolume_file(iter(OrientableCuspedCensus),
                         PSEUDO_VOLUME_FILENAME,
                         max_secs = 160,
                         sln = n,
                         engine = None,
                         retrieve = True)

# Read in pseudo-volumes
pseudo_volumes = read_volumedata_csv(PSEUDO_VOLUME_FILENAME)
print 'Read pseudovolume file'
pseudo_volumes.clean(n=n, maxsfdegree=maxdegree)
print 'Cleaned pseudovolume file'

# Read in span data
spans = SpanData(read_spans('spans.csv'))
print 'Read span file'
spans.fit(pseudo_volumes, n=n)
print 'Fitted pseudovolumes to spans'

# Read in all the other data
dataset = quick_read_csv('all_volumes.csv')
print 'Read in all volumes'

# spans.write_to_csv('something', dataset)
spans.write_nice_fits('linear_combinations_for_n='+str(n)+'.csv')
print 'Wrote nice fits'
spans.write_failures('linear_combination_fit_failures_for_n='+str(n)+'.csv')
