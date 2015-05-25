from cypari import *
import itertools
import re
import time
import tempfile
import os.path

from VolumeUtilities import *
from VolumeFinder import _kickoff_snap, _send_cmd, _drain
from VolumeFinder import RE_INV_TRACE_FIELD_NOT_FOUND, RE_FUNC_REQ_GROUP, RE_ERROR

_RE_NO_BOREL = re.compile('.*Borel Regulator: not computed.*', re.DOTALL)
_RE_BOREL_REGULATOR = re.compile('.*Borel Regulator: \[([-0-9.,]*)\].*', re.DOTALL)
_RE_BAD_GENERATOR = re.compile('.*Non existent generator.*', re.DOTALL)

def find_span(elts, n, cutoff = 4096):
    """
    Parameters: elts is a list [a1, a2, a3, ... ], where each ai is a
    (volume (as pari object), name), and n is the upper limit of the
    dimension of the span to be found.

    Returns: returns (s, sn, r, bm, bd), where

    s is determined to be the optimal span that is a subset of elts for
    an n-dimensional lattice,

    sn is the corresponding names for s,

    r is the ratio of det(S) / d, with S the matrix formed out of s, and
    d the gcd of all possible such determinants (it follows that s is a
    basis for elts iff r is 1),

    Optionally, the parameter cutoff may be passed (defaulting to
    4096). This controls the upper bound on absolute values of
    coefficients resulting from pari's lindep that are considered
    relevant.  For example [sqrt(2), 0.00001] gives a linear dependence
    of [-160, 22627417]. In practice, this should probably be ignored.

    Note that if a span that would have greater dimension than n would
    be generated, this method will raise a ValueError. This indicates
    that the input must be invalid.
    """

    # Remove integer multiples, creates something like [2, 5, 7, 9, 11, ... ]
    # This could be done in the next step, but is singled out as a special case
    # of the general operation.
    base_elts = sorted(elts, key = lambda a: -abs(float(a[0])))
    base_elts = [(e[0].replace('-',''), e[1]) for e in base_elts]   # take absolute value
    independent_elts = []
    while base_elts:
        cur_min_elt = base_elts.pop()
        independent_elts.append(cur_min_elt)
        base_elts = [ a for a in base_elts if not _is_int(float(a[0])/float(cur_min_elt[0])) ]

    # Create a rational basis, but be careful of [ 0.5, 0.75, sqrt(2) ], which
    # might select first two as e1 and e2, yet gives pari no representation of
    # sqrt(2) by the other two
    next_pure_basis_element_position = 0 # shall create [1, 0, ..., 0] as next ei
    rational_vectors = None # this will be a list of tuples, (vol, vec, name), with vec in Q^n
    rational_basis = None # shall be { vol where vec is some ei }
    for i in range(0, len(independent_elts)):
        vol, name = independent_elts[i] # gen.pari(independent_elts[i])
        if not rational_basis:
            next_pure_basis_element_position = 1
            rational_basis = [vol]
            rational_vectors = [(vol,gen.pari([1]).mattranspose(),name)]
        else:
            test_vec = rational_basis[:]
            test_vec.insert(0, vol)
            lindep_results = gen.pari(test_vec).lindep(LINDEP_PRECISION)
            if lindep_results[0] == 1:
                # This volume may be expressed as a linear, integral
                # combination of others. It is worthless.
                pass
            elif lindep_results[0] == 0 or max(abs(lindep_results)) > cutoff:
                # This volume is considered by pari to be harder to fit into
                # the current set of results than allowable. Treat it as a new
                # basis element.
                if next_pure_basis_element_position >= n:
                    raise ValueError('The current basis elements are ' + str(rational_basis) +
                                     ', and attempting to fit volume ' + str(vol) +
                                     ' into them determined it should be classed as a new basis element, but there are already ' +
                                     str(n) + ' such!')
                vec = [0] * (next_pure_basis_element_position + 1)
                vec[next_pure_basis_element_position] = 1
                next_pure_basis_element_position += 1
                rational_vectors = [ (c_p,c_v.concat(gen.pari([0])),c_n) for c_p,c_v,c_n in rational_vectors ]
                rational_basis.append(vol)
                rational_vectors.append((vol,gen.pari(vec).mattranspose(),name))
            else:
                # This volume is a rational combination of the selected basis
                # elements. Use it.
                vec = - lindep_results[1:] / lindep_results[0]
                while len(vec) < next_pure_basis_element_position:
                    vec = vec.concat(gen.pari([0]))
                rational_vectors.append((vol,vec.mattranspose(),name))

    # Now rational_vectors is a python list containing pari row vectors representing each elt. Get dets and such
    determined_n = next_pure_basis_element_position # if there is less data than ncp, keep going
    det_gcd = None
    min_det = None
    min_det_elts = None
    min_det_names = None
    for mat_tuple in itertools.combinations(rational_vectors, determined_n):
        mat_elts = [vec for vol,vec,name in list(mat_tuple)]
        names = [name for vol,vec,name in list(mat_tuple)]
        mat = mat_elts[0]
        for v in mat_elts[1:]:
            mat = mat.concat(v)
        cur_det = mat.matdet().abs()

        if cur_det == 0:
            continue

        if det_gcd is None:
            det_gcd = cur_det
        else:
            det_gcd = det_gcd.gcd(cur_det)

        if (min_det is None) or (cur_det.numerator()*min_det.denominator() < min_det.numerator()*cur_det.denominator()):
            min_det = cur_det
            min_det_elts = [vol for vol,vec,name in list(mat_tuple)]
            min_det_names = names

    return min_det_elts, min_det / det_gcd, min_det_names

def find_borel_matrix(manifold_names, shape_field_degree):
    """
    Given manifold_names, a list of strings each the name of a
    manifold (i.e., they may be passed as a constructor to SnapPy's
    Manifold object) returns (borel_matrix, borel_matrix_determinant).

    borel_matrix is a list of the Borel regulators (themselves lists) for
    each manifold in sn (not a Pari matrix object)

    borel_matrix_determinant is the determinant (a pari object) of bm if
    it can be formed into a square matrix (all Borel regulators have the
    same dimension, and this dimension is the same as that of
    manifold_names), or None if it cannot.

    The parameter shape_field_degree is passed directly to
    snap as by the command `set degree'. As it increases, the
    computations become slower, but fewer exceptions are raised.

    Note: This method currently relies on the external program `snap` to
    compute borel regulators. If errors arise, try adjusting the options
    for snap found in VolumeUtilities.
    """
    borel_reg_strings = list()

    next_borel_reg = _borel_regulator(manifold_names[0], shape_field_degree)
    borel_reg_mat = gen.pari(next_borel_reg).mattranspose()
    borel_reg_strings.append(next_borel_reg)

    for next_name in manifold_names[1:]:
        next_borel_reg = _borel_regulator(next_name, shape_field_degree)
        borel_reg_mat = borel_reg_mat.concat(gen.pari(next_borel_reg).mattranspose())
        borel_reg_strings.append(next_borel_reg)

    borel_determinant = None
    try:
        if len(borel_reg_mat) == len(borel_reg_mat[0]):
            borel_determinant = borel_reg_mat.matdet()
    except:
        pass

    return borel_reg_strings, borel_determinant

def _borel_regulator(manifold, shape_field_degree, temp_dir = None):
    """
    Given a manifold (either its name as a string, or a Manifold
    object itself), return, as a list, the volumes (in string form)
    that make up that manifold's Borel regulator. The parameter
    shape_field_degree is passed directly to snap as by the command
    `set degree'. As it increases, the computations become slower, but
    fewer exceptions are raised.

    Ex: _borel_regulator('m349') returns ['2.721625...', '4.666479...'].

    Note: By default, for each invocation of this method a temporary
    directory is created. To avoid clutter, the optional argument
    temp_dir can be passed instead. This may introduce race conditions
    if other threads are using the same temp_dir for this method, but
    if such a temp_dir is known to be used by only one thread, this
    method will remain threadsafe.
    """
    if type(manifold) is str:
        input_manifold = Manifold(manifold)
    else:
        input_manifold = manifold

    if temp_dir:
        trig_file_dir = temp_dir
    else:
        trig_file_dir = tempfile.mkdtemp()

    input_manifold.save(os.path.join(trig_file_dir, 'tmp.trig'))
    snap_process = _kickoff_snap(trig_file_dir)
    try:
        _send_cmd(snap_process, 'set precision 25')
        _send_cmd(snap_process, 'set degree ' + str(int(shape_field_degree)))
        _send_cmd(snap_process, 'read file tmp.trig')
        _send_cmd(snap_process, 'co inv co sh print ari')
        snap_output = ''
        wait_count = 0
        while True:
            snap_output += _drain(snap_process.stdout)
            borel_match = _RE_BOREL_REGULATOR.match(snap_output)
            if borel_match is not None:
                return re.split(',', borel_match.group(1))
            elif any([RE_INV_TRACE_FIELD_NOT_FOUND.match(snap_output),
                        RE_FUNC_REQ_GROUP.match(snap_output),
                        RE_ERROR.match(snap_output),
                        _RE_BAD_GENERATOR.match(snap_output),
                        _RE_NO_BOREL.match(snap_output)]):
                raise Exception('Manifold ' + str(input_manifold) + ' has incalculable Borel regulator')
            else:
                time.sleep(0.25)
                wait_count += 1
                if wait_count > 50:
                    raise Exception('Manifold ' + str(input_manifold) + ' appears to have hung snap')
    finally:
        snap_process.terminate()

def _is_int(r):
    s = r % 1
    return s < EPSILON or s + EPSILON > 1
