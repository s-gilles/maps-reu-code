from snappy import *
from math import log

# Global constants:
PARI_PRECISION = 500
LINDEP_PRECISION = 16
EPSILON = 1e-12
MAX_ITF = 8
MAX_COEFF = 4096
MAX_LDP_TRIES = 4
MAX_LINDEP_TRIES = 50

# This is designed to match SnapPy's solution_type(enum = True) output
SOL_TYPE_STRINGS = ['not_attempted', 'geometric', 'nongeometric', 'flat', 'degenerate', 'unrecognized', 'none_found']


# *** ATTENTION USER ***
#
# If you want to dynamicallly change the above constants for all
# sessions, or just want to keep the default values so you can revert,
# do so here! That way the values will both be used and tested in the
# warning system below

# User definitions:

# Enviroment setup code:
pari.set_real_precision(PARI_PRECISION)

# Test code; makes sure the constants are are sanely set:
# PARI & LINDEP PRECISION
if .9*PARI_PRECISION <= LINDEP_PRECISION or PARI_PRECISION - 3 <= LINDEP_PRECISION:
    print 'WARNING: You set PARI to use '+str(PARI_PRECISION)+' places by default, with lindep calls at '+str(LINDEP_PRECISION)+' places;'
    print 'This will probably read to rounding errors messing things up when lindep is used.'
    print 'You should probably make sure LINDEP_PRECISION < both .9*PARI_PRECISION and PARI_PRECISION - 3 to avoid this.'
# EPSILON (vs. PARI_PRECISION)
if EPSILON <= 0:
    print 'WARNING: You set EPSILON to '+str(EPSILON)+', but it must be positive.'
    print 'Try setting EPSILON=abs(EPSILON)'
if EPSILON > .01:
    print 'WARNING: You set EPSILON to '+str(EPSILON)+', which is really big.'
    print 'PARI is capable of computing to hundreds of places of PRECISION, currently '+str(PARI_PRECISION)+', and at worst will use over 10.'
    print 'You should set EPSILON to something smaller.'
if log(EPSILON) >= .9*log(PARI_PRECISION):
    print 'WARNING: You set EPSILON to '+str(EPSILON)+', which is small compared to PARI\'s PRECISION limit of '+str(float(10)**(-PARI_PRECISION))+'.'
    print 'Maybe you should set EPSILON a little bigger, so rounding errors from PARI don\'t get through.'
