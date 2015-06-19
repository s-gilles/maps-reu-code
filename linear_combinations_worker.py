#!/usr/bin/python2

from __future__ import print_function

import argparse
import errno
import fractions
import csv
import sys
import os

from snappy import *

# This is the `disposable' client.  Run it with args as specified (run
# without args for options) stdin, it accepts a filename. It then
# expects to receive manifold names, one afterthe other, on stdin. It
# will try to find the spanning of that manifold's dimensions
# (geometric and non-geometric) and output them to the filename
# (headers should be handled correctly).

LCC_READY_TO_RECEIVE_DATA_STRINGS = 'LCC_I_am_ready_for_data'
LCC_MAKING_PROGRESS = 'LCC_I_am_making_progress'
LCC_WANT_MORE_DATA = 'LCC_ready_for_more_data'
LCC_ENCOUNTERED_ERROR = 'LCC_encountered_an_error'
LCC_REFUSE_TO_CLOBBER = 'LCC_output_file_exists_and_is_not_what_was_expected'
LCC_CANT_DEAL_WITH_OUTPUT_FILE = 'LCC_output_file_cannot_be_gracefully_dealt_with'
LCC_HEADER = 'Manifold;TraceField;TraceFieldDegree;Subfield;SubfieldDegree;Root;Volume;LinearCombination'

MAX_COEFF = 4096
MAX_LINDEP_TRIES = 50

class Span:

    polynomial_string = None
    """The string defining the field of this span. Example 'x^7-3*x^6+5*x^5-5*x^4+2*x^3-2*x+1'"""

    number_of_complex_places = None
    """The number of complex places of polynomial_string. Example 2"""

    root = None
    """Each span represents a specific complex root of the polynomial with positive imaginary coefficient. Example '0.123338470+1.00512926*I'"""

    span_dimension = None
    """How many elements are in the volume span. Example: 1"""

    volume_span = None
    """A list, with span_dimension elements, which defines a linearly
    independent subset of the volumes detected for this polynomial and
    root.

    Example: ['7.53260450387045224483849916956535836']

    """

    manifold_span = None
    """A list, with span_dimension elements, such that each manifold has
    the geometric volume of the corresponding entry in volume_span.

    Example: ['10^3_17(-1,1)(4,1)(3,1)']

    """

    fit_ratio = None
    """An integer representing how well the manifold_span matched observed
    data. If, for example, the lattice of detected volumes contains
    9*x and 10*x, yet the smallest volume found is 5*x, the (if the
    lattice is really generated by geometric volumes of manifolds)
    there must be a manifold with volume x, yet such could not be
    obtained - the best is a manifold with volume 5*x. So the
    fit_ratio would be 5.

    Example: 1

    """


    def __init__(self, _polynomial, _degree, _number_of_complex_places, _root,
                 _span_dimension, _volume_span, _manifold_span,
                 _fit_ratio, _borel_regulator_matrix, _borel_regulator_det):
        """Constructor"""
        self.polynomial = _polynomial
        self.number_of_complex_places = int(_number_of_complex_places)
        self.root = _root
        self.span_dimension = int(_span_dimension)
        self.volume_span = eval(_volume_span)
        self.manifold_span = eval(_manifold_span)
        self.fit_ratio = int(_fit_ratio)

        # For speed, don't check anything here, although we could
        # assert that len(volume_span) == span_dimension, that
        # everything in manifold_span is a valid manifold, etc. etc.


all_spans = []
"""A list of every span read in.  Dunno why"""

span_dict = dict()
"""A lookup of spans by polynomial"""

output_filename = None
"""Where to write out data"""

n_in_SLnC = 2
"""Which SL(n,C) to consider representations into"""

lindep_precision = 16
"""How much precision to ask Pari to use"""

epsilon_string = '1E-500'
"""Very small number for use in volume culling"""

nfsubfields_cache = dict()
"""Cache for the costly job of computing subfields"""

def get_nfsubfields(polynomial):
    if polynomial in nfsubfields_cache:
        return nfsubfields_cache[polynomial]
    done = False
    nfsubfields_output = []
    polredabs_output = None
    while True:
        try:
            if not polredabs_output:
                try:
                    polredabs_output = pari(polynomial).polredabs()
                except:
                    polredabs_output = pari(polynomial)
                sys.stdout.flush()
                nfsubfields_output = polredabs_output.nfsubfields()[1:]
            break
        except:
            pari.allocatemem()

    nfsubfields_output = [ str(x[0]) for x in nfsubfields_output
                           if len(str(x[0])) > 2]
    sys.stdout.flush()
    reduced = list()
    for subfield in nfsubfields_output:
        while True:
            try:
                reduced.append(str(pari(subfield).polredabs()))
                sys.stdout.flush()
                break
            except:
                pari.allocatemem()
    reduced = [ s.replace(' ', '') for s in reduced ]
    nfsubfields_cache[polynomial] = reduced
    return reduced

def get_pari_lindep(str_vols,
                    maxcoeff = MAX_COEFF,
                    max_tries = MAX_LINDEP_TRIES):
    """Given str_volumes, a list of volumes in string form, returns the
    dependancy found (if any) as a list of integers if all coefficents
    are <= maxcoeff or maxcoeff is nonpositive; otherwise, it returns
    []

    """
    vols = list(str_vols)   # in case someone sent some other type collection
    vec = None

    num_tries = 0
    while not vec and num_tries < max_tries:
        vec = str(pari(str(vols).replace("\'",'')).lindep(lindep_precision))[1:-2].replace(' ','').split(',')
        num_tries += 1

    if not vec or vec == ['']: # no input
        print(LCC_ENCOUNTERED_ERROR)
        sys.stdout.flush()
        return list()

    o = [int(v) for v in vec]
    if len(o) >= 1 and o[-1] < 0:
        o = [-1 * element for element in o]
    if maxcoeff > 0:
        for x in o:
            if abs(x) > maxcoeff:
                return list()
                break
    return o

def handle_manifold(a_string):
    m = Manifold(a_string)
    pv = m.ptolemy_variety(n_in_SLnC, 'all')
    decompositions = None
    try:
        decompositions = pv.retrieve_decomposition()
    except:
        print(LCC_ENCOUNTERED_ERROR)
        sys.stdout.flush()
        return
    subfields = None
    spans_for_polynomial = None

    for decomposition in decompositions:
        obstruction_class = -1
        spans_for_polynomial = None
        for s in decomposition:
            obstruction_class = obstruction_class + 1
            polynomial = str(s.number_field()).replace(' ', '')
            subfields = None
            spans_for_polynomial = None
            volumes = s.solutions(numerical = True).volume_numerical()
            if not volumes:
                volumes = []

            # First, prune tiny volumes
            volumes = [ str(v) for v in volumes if pari('abs('+str(v)+')>'+epsilon_string) ]

            # Second, ensure that there are not more than 4 = 8/2
            # volumes, as this would indicate a field we aren't
            # prepared to deal with.
            distinct_volumes = list()
            for v in volumes:
                matches = [ u for u in distinct_volumes if pari('abs(abs('+str(v)+')-abs('+str(u)+'))<='+epsilon_string) ]
                if not matches:
                    distinct_volumes.append(v)
            if len(distinct_volumes) > 4:
                continue

            # Now we need to get the possible subfields. This is the
            # most intensive part.
            if not subfields:
                subfields = get_nfsubfields(polynomial)
                print(LCC_MAKING_PROGRESS)
                sys.stdout.flush()

            for subfield in subfields:
                # For every root of the polynomial, try and fit the data
                # into the span
                spans_for_polynomial = span_dict.setdefault(subfield, list())

                for volume in distinct_volumes:
                    smallest_fit = -1
                    best_span = None
                    ldp = []

                    # Check all spans to get the one with the best fit
                    # ratio for this volume
                    for span in spans_for_polynomial:
                        ldp = get_pari_lindep(span.volume_span + [ volume ])
                        if not ldp or ldp[-1] == 0:
                            continue
                        if smallest_fit < 0 or ldp[-1] < smallest_fit:
                            smallest_fit = ldp[-1]
                            best_span = span
                        if ldp[-1] == 1:
                            break

                    if not ldp or ldp[-1] == 0:
                        continue

                    # Note we lie a little here because we want to
                    # make sure that we aren't killed while writing
                    # the output
                    print(LCC_MAKING_PROGRESS)
                    sys.stdout.flush()

                    linear_comb_str = ''
                    for n in xrange(len(ldp)-1):
                        if ldp[n] == 0:
                            continue
                        if linear_comb_str != '' and ldp[n] < 0:
                            linear_comb_str += ' + '
                        linear_comb_str += str(fractions.Fraction(-1 * ldp[n], ldp[-1])).replace('-', ' - ')
                        linear_comb_str += '*'
                        linear_comb_str += best_span.manifold_span[n]

                    with open(output_filename, 'ab') as output_file:
                        csvwriter = csv.writer(output_file,
                                               delimiter = ';',
                                               quoting = csv.QUOTE_ALL,
                                               lineterminator = os.linesep)
                        csvwriter.writerow([str(m),
                                            str(polynomial),
                                            str(pari(polynomial).poldegree()),
                                            str(subfield),
                                            str(pari(subfield).poldegree()),
                                            str(best_span.root),
                                            str(volume),
                                            linear_comb_str.strip()])

if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('--spanfile', '-s',
                            default = 'spans.csv',
                            help = 'Path of master list of spans for reading')
    arg_parser.add_argument('--output', '-o',
                            default = 'linear_combinations.csv',
                            help = 'Output path (appending handled well)')
    arg_parser.add_argument('-n',
                            default = '2',
                            metavar = 'n',
                            help = 'For considering SL(n,C)')
    arg_parser.add_argument('--precision', '-p',
                            type = int,
                            default = 500,
                            help = 'Precision to set PARI to')
    arg_parser.add_argument('--lindep-precision', '-l',
                            type = int,
                            default = 16,
                            help = 'Precision to use when calling lindep()')
    args = arg_parser.parse_args()
    output_filename = args.output
    n_in_SLnC = args.n
    if args.precision < 20:
        args.precision = 20
    pari.set_real_precision(args.precision)
    epsilon_string = '1E-'+str(args.precision-3)
    lindep_precision = args.lindep_precision

    # Figure out how to handle the output file - should we append the
    # header or not?
    has_header = False
    must_abort = False

    try:
        with open(output_filename, 'rb') as output_file:
            firstline = output_file.readline().rstrip()
            if firstline:
                if firstline == LCC_HEADER:
                    has_header = True
                else:
                    print(LCC_REFUSE_TO_CLOBBER)
                    sys.stdout.flush()
                    must_abort = True
    except IOError as e:
        if e.errno == errno.ENOENT:
            pass
        else:
            print(LCC_CANT_DEAL_WITH_OUTPUT_FILE)
            sys.stdout.flush()
            must_abort = True
    except:
        print(LCC_CANT_DEAL_WITH_OUTPUT_FILE)
        sys.stdout.flush()
        must_abort = True

    if must_abort:
        sys.exit(1)

    if not has_header:
        with open(output_filename, 'a') as output_file:
            output_file.write(LCC_HEADER)
            output_file.write('\n')

    # Now the output file has the header. Note that a handle is not
    # kept because this process could be killed at any time, so as
    # close to atomicity as possible is desired.

    with open(args.spanfile, 'rb') as span_file:
        reader = csv.reader(span_file, delimiter=';', quoting = csv.QUOTE_ALL)
        header = reader.next()
        for row in reader:
            new_span = Span(*row)
            all_spans.append(Span(*row))

            spans_for_poly = span_dict.setdefault(new_span.polynomial, list())
            spans_for_poly.append(new_span)


    print(LCC_READY_TO_RECEIVE_DATA_STRINGS)
    sys.stdout.flush()

    while True:
        handle_manifold(sys.stdin.readline())
        print(LCC_WANT_MORE_DATA)
        sys.stdout.flush()
