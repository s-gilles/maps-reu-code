from cypari import *
from itertools import product
from operator import __mul__

def _matker(elts):
    """
    Given an array of elements [a1, a2, ... an], where each ai is a pari
    polymod object (assumed to be over the same polynomial, return a
    linear dependency between the lifts of each ai in the form of a
    matrix kernel. This will be returned as an empty matrix if none such
    can be found (which evaluates to False when converted to a boolean).

    For example, if a, b, and c represent x^2, x + 1, and x^2 - x - 1,
    _matker([a, b, c]) might return [1, -1, -1] as the result. If the
    matrix produced by taking the vectors of each [lifted] polynomial
    admits multiple kernels, each will be returned as an entry in the
    matrix.
    """
    n = 1 + max(p.lift().poldegree() for p in elts)
    elt_vecs = list()
    for e in elts:
        v = e.lift().Vec()
        padding_length = n - len(v)
        v = [0] * padding_length + list(v)
        elt_vecs.append(gen.pari(v).mattranspose())

    m = gen.pari(elt_vecs[0])
    for next_elt in elt_vecs[1:]:
        m = m.concat(next_elt)

    return m.matker()

def _all_factors_of(n):
    """
    Given an integer n, return, in list form, the set {d : 1 <= d <= n
    and d divides n }. This may be used to simulate PARI's divisors()
    function, which is not exposed at the time of this writing.
    """
    frange = [[f[0] ** e for e in range(0,f[1]+1)] for f in gen.pari(n).factor().mattranspose()]
    all_factors = [reduce(__mul__, list(i), 1) for i in product(*frange)]
    return [f for f in all_factors if f < n ]

def _degree_over_Q(x):
    """
    The minimal n such that there is a polynomial p of degree n where
    p(x) = 0
    """
    dep_list = [gen.pari(1).Mod(x.mod())]
    n = 0
    while True:
        n += 1
        dep_list.insert(0, x ** n)
        if bool(_matker(dep_list)):
            return n

def _degree_over_Q_adjoined(x, a):
    """
    The minimal n such that there is a polynomial of degree n in Q(a)
    where p(x) = 0
    """
    dep_list = [a, gen.pari(1).Mod(x.mod())]
    n = 0
    while True:
        n += 1
        dep_list.insert(0, x ** n)
        dep_list.insert(0, a * (x ** n))
        if bool(_matker(dep_list)):
            return n

class _SmallIntegerIterator:
    """
    An iterator that returns 1, -1, 2, -2, 3, -3, ...
    """
    def __init__(self, max_val = -1):
        self.next_int = 1
        self.max_val = max_val

    def __iter__(self):
        return self

    def next(self):
        last_int = self.next_int

        if self.max_val > 0 and last_int > self.max_val:
            raise StopIteration

        if self.next_int > 0:
            self.next_int = -self.next_int
        else:
            self.next_int = -self.next_int + 1
        return last_int

def reduce_elements(alpha,
                    elts,
                    max_coefficient = -1):
    """
    Given

    - alpha, a pari polymod (e.g. Mod(x, x^2 - x - 1))

    - elts, a list [b1, b2, ... bn], where each bi is a pari polymod
    element that is itself a polynomial of alpha

    computes a gamma such that Q(gamma) is Q(b1, b2, ... bn), and
    returns

    - gamma such that Q(gamma) = Q(b1, b2, ... bn)

    - q, the minimum polynomial of gamma

    - a list [k1, k2, ... kn] such that gamma = Sum (ki * bi)

    - a list [c1, c2, ... cn], where each ci is a Mod(q, r), such that r
    is the minimum polynomial of gamma, and q is a polynomial such that
    q(gamma) = bi: in other words the list elts modified to be over
    Q(gamma), not Q(alpha).

    If gamma cannot be computed (e.g., if one of elts is not over
    alpha), an error is raised.

    The optional parameter max_coefficient restricts the absolute value
    of each ki. If max_coefficient is specified, but is too low for the
    algorithm to append any ki * bi for a bi which has been determined
    useful, a ValueError will be raised. This will, however, prevent
    possibly unbounded execution time in the case of malformed input.
    """

    # Calculate [Q(alpha) : Q]. This can be used to break out of the
    # beta loop if it is determined that Q(gamma) = Q(alpha)
    alpha_degree = 0
    dep_list = [gen.pari(1).Mod(alpha.mod())]
    for alpha_degree in range(1,alpha.mod().poldegree() + 2):
        dep_list.insert(0,alpha ** alpha_degree)
        if bool(_matker(gen.pari(dep_list))):
            break

    # Initialize the coefficients [k1, k2, ... kn] for gamma's formal
    # sum
    gamma_coefficients = [0] * len(elts)

    gamma = elts[0]
    gamma_degree = _degree_over_Q(gamma)

    i = 0
    gamma_coefficients[i] = 1

    # Adjoin each beta element, one at a time, increasing the degree
    # of the extension as much as possible. If any adjoinment raises
    # the degree of the extension to alpha's degree, the adjoinment
    # must be complete, so the loop is short-circuited
    for beta in iter(elts[1:]):
        i += 1

        # If the loop was exhausted and gamma was updated, check if
        # the degree of gamma is now the degree of alpha. If so, no
        # futher beta can increase the extension's degree, so return.
        if gamma_degree >= alpha_degree:
            break

        # [ Q(b, g) : Q ] = [ Q(b, g) : Q(g) ] * [ Q(g) : Q ]
        m = _degree_over_Q_adjoined(beta, gamma)
        if m == 1:
            continue

        n = m * gamma_degree

        # Build { d : d | n, 0 < d < n } for use in testing degrees
        all_factors = _all_factors_of(m)

        # Pick a smallish r
        r_choice_works = False
        for r in _SmallIntegerIterator(max_coefficient):
            proposed_gamma = gamma + beta * r
            r_choice_works = True

            # This next part could be a check of
            # _degree_over_Q_adjoined(proposed_gamma, gamma) against m
            # = _degree_over_Q_adjoined(beta, gamma), but since it
            # must be a divisor of m, this is marginally faster.

            powers_of_gamma = list()
            powers_of_gamma.append(gen.pari(1).Mod(alpha.mod()))
            powers_of_gamma.append(gamma)
            for d in range(1, m+1):
                powers_of_gamma.append(proposed_gamma ** d)
                powers_of_gamma.append(gamma * (proposed_gamma ** d))

            # Make sure [ Q(proposed_gamma, gamma) : Q(gamma) ] <= m
            if not bool(_matker(powers_of_gamma)):
                r_choice_works = False

            # Make sure [ Q(proposed_gamma, gamma) : Q(gamma) ] >= n
            for f in all_factors:
                if r_choice_works:
                    if bool(_matker(powers_of_gamma[0:2 * (f + 1)])):
                        r_choice_works = False

            # By this point, it should be that
            # r_choice_works = (_degree_over_Q_adjoined(proposed_gamma, gamma) == m)
            # which should be equivalent to
            # r_choice_works = (_degree_over_Q(proposed_gamma) == n)

            # If the r works, simply update gamma and leave
            if r_choice_works:
                gamma = proposed_gamma
                gamma_degree = n
                gamma_coefficients[i] = r
                break

        # If the loop was exhausted with no r found, leave - the data
        # was strange
        if not r_choice_works:
            return ValueError

    # Now that gamma is computed, each beta in elts has to be
    # recomputed over gamma instead of alpha. First, minpoly(gamma) is
    # computed, then a polynomial for each beta is built up out of it.
    gamma_minpoly = gen.pari('minpoly(' + str(gamma) + ')')
    x_over_gamma = gen.pari('x').Mod(gamma_minpoly)
    beta_primes = list()
    for beta in elts:
        dep_list = [beta, gen.pari(1).Mod(alpha.mod())]
        d = 0
        dep_results = None
        while True:
            d += 1
            dep_list.append(gamma ** d)
            dep_results = _matker(dep_list)
            if bool(dep_results):
                beta_prime = gen.pari('0').Mod(gamma_minpoly)
                break

        d = 0
        for k in dep_results[0][1:]:
            beta_prime += k * (x_over_gamma ** d)
            d += 1
        beta_primes.append(beta_prime * gen.pari(-1) / gen.pari(dep_results[0][0]))

    return gamma, gamma_minpoly, gamma_coefficients, beta_primes
