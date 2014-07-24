from cypari import *
import itertools

from VolumeUtilities import *

def find_span(elts, n, cutoff = 4096):
    """Parameters: elts is a list [a1, a2, a3, ... ], where each ai is a (volume (as pari object), name)

Returns: returns (s, sn, r), where s is determined to be the optimal span that is a subset of elts for an n-dimensional lattice, sn is the corresponding names for s, and r is the ratio of det(S) / d, with S the matrix formed out of s, and d the gcd of all possible such determinants.  It follows that s is a basis for elts iff r is 1.

Optionally, the parameter cutoff may be passed (defaulting to 4096). This controls the upper bound on absolute values of coefficients resulting from pari's lindep that are considered relevant.  For example [sqrt(2), 0.00001] gives a linear dependence of [-160, 22627417]. In practice, this should probably be ignored."""

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

def _is_int(r):
    s = r % 1
    return s < EPSILON or s + EPSILON > 1
