"""Exact local branch recovery on the published Paparazzi trefoil.

Geometry: Yao et al., IEEE T-RO 37(4), 2021, Section VI-D,
doi:10.1109/TRO.2020.3043690. An endpoint cone and rate-duration guard
certify coverage by both local charts. Other regions return OUTSIDE_COVERAGE.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import ceil, floor

from monotone_history import preimage
from straight_history import intersect, propagate, witness

CHARTS = {'negative': (F(-11, 10), F(-9, 10)),
          'positive': (F(9, 10), F(11, 10))}
RATES = (F(27, 1000), F(33, 1000))
LOCAL_PRIOR = 'both-local-crossing-charts'
SOURCE = 'https://arxiv.org/pdf/2012.01826v3'


def curve(q):
    """Rational half-angle substitution q=tan(u/2) into the source equation."""
    q = F(q)
    d = 1+q*q
    c, s = (1-q*q)/d, 2*q/d
    radius = 160+80*(4*c*c*c-3*c)
    return radius*(c*c-s*s)+79, radius*2*s*c-F(681, 10)


def _multiply(left, right):
    result = [0]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] += a*b
    return tuple(result)


# Numerators in increasing powers, over (1+q^2)^5, before translation.
NUMERATORS = tuple(tuple(80*c for c in _multiply((3, 0, -9, 0, 21, 0, 1), p))
                   for p in ((1, 0, -6, 0, 1), (0, 4, 0, -4)))


def _horner(coefficients, q):
    value = F(0)
    for coefficient in reversed(coefficients):
        value = value*q+coefficient
    return value


def expanded_curve(q):
    """Independent expanded equation used by the scorer/witness checker."""
    q = F(q)
    d = (1+q*q)**5
    return (_horner(NUMERATORS[0], q)/d+79,
            _horner(NUMERATORS[1], q)/d-F(681, 10))


def _derivative_numerator(coefficients):
    # d/dq [N(q)/(1+q^2)^5] = [N'(q)(1+q^2)-10qN(q)]/(1+q^2)^6.
    derivative = tuple(k*coefficients[k] for k in range(1, len(coefficients)))
    result = list(_multiply(derivative, (1, 0, 1)))
    result.extend([0]*(len(coefficients)+1-len(result)))
    for k, coefficient in enumerate(coefficients):
        result[k+1] -= 10*coefficient
    return tuple(result)


DERIVATIVES = tuple(_derivative_numerator(p) for p in NUMERATORS)


def _interval_product(a, b):
    products = tuple(x*y for x in a for y in b)
    return min(products), max(products)


def _interval_horner(coefficients, interval):
    result = (F(0), F(0))
    for coefficient in reversed(coefficients):
        low, high = _interval_product(result, interval)
        result = low+coefficient, high+coefficient
    return result


def derivative_enclosure(domain, axis):
    """Rational outward enclosure on one sign-definite parameter interval."""
    low, high = domain
    if low > high or low <= 0 <= high or axis not in (0, 1):
        raise ValueError('ordered sign-definite interval and coordinate required')
    squares = sorted((low*low, high*high))
    denominator = ((1+squares[0])**6, (1+squares[1])**6)
    numerator = _interval_horner(DERIVATIVES[axis], domain)
    return _interval_product(numerator, (1/denominator[1], 1/denominator[0]))


@lru_cache(maxsize=1)
def geometry_certificate():
    """Prove monotonicity by 256 rational cells per chart, without sampling.

    Interval Horner includes every point of each closed cell. Positive
    denominators and their reciprocals are enclosed exactly. Integer rounding
    is outward; the resulting bounds can be checked again by the verifier.
    """
    subdivisions, bounds = 256, {}
    for name, (low, high) in CHARTS.items():
        axis_bounds = []
        for axis in (0, 1):
            enclosures = [derivative_enclosure(
                (low+(high-low)*k/subdivisions,
                 low+(high-low)*(k+1)/subdivisions), axis)
                for k in range(subdivisions)]
            bound = (F(floor(min(v[0] for v in enclosures))),
                     F(ceil(max(v[1] for v in enclosures))))
            if bound[0] <= 0 <= bound[1]:
                raise ArithmeticError('coordinate monotonicity not certified')
            axis_bounds.append(bound)
        bounds[name] = tuple(axis_bounds)
    return dict(subdivisions=subdivisions, derivative_bounds=bounds,
                method='rational interval Horner on a complete closed-cell partition')


def verify_geometry_certificate(certificate):
    count = certificate['subdivisions']
    if type(count) is not int or count < 1:
        return False
    for name, (low, high) in CHARTS.items():
        for axis in (0, 1):
            claimed = tuple(map(F, certificate['derivative_bounds'][name][axis]))
            if claimed[0] <= 0 <= claimed[1] or claimed[0] > claimed[1]:
                return False
            for k in range(count):
                bound = derivative_enclosure((low+(high-low)*k/count,
                                               low+(high-low)*(k+1)/count), axis)
                if not claimed[0] <= bound[0] <= bound[1] <= claimed[1]:
                    return False
    return True


def parse(raw):
    times = tuple(map(F, raw['times']))
    positions = tuple(tuple(map(F, row)) for row in raw['positions'])
    errors = tuple(map(F, raw['linf_errors']))
    if (not times or len(times) != len(positions) or len(times) != len(errors)
            or any(len(row) != 2 for row in positions) or any(e < 0 for e in errors)
            or any(a >= b for a, b in zip(times, times[1:]))):
        raise ValueError('nonempty equal-length 2D histories with strict timestamps required')
    return times, positions, errors


def verify_witness(values, raw, domain):
    """Direct expanded-coordinate/rate check; does not use curve or inverses."""
    times, positions, errors = parse(raw)
    values = tuple(map(F, values))
    return (len(values) == len(times)
            and all(domain[0] <= q <= domain[1] for q in values)
            and all(abs(x-y) <= error for q, point, error in zip(values, positions, errors)
                    for x, y in zip(expanded_curve(q), point))
            and all(RATES[0]*(v-u) <= q-p <= RATES[1]*(v-u)
                    for u, v, p, q in zip(times, times[1:], values, values[1:])))


def coverage_guard(raw):
    """Sufficient full-history chart coverage inferred from measurements.

    The source polar radius is positive. Endpoint error boxes inside the cone
    X<0, |Y/(-X)|<1/6 put every compatible endpoint inside the two charts.
    Their boundary cone ratios exceed 1/6. The rate cap excludes crossing
    the chart gap; positive timing then confines intermediate parameters.
    An optional legacy 'prior' field neither grants nor bypasses this guard.
    """
    times, positions, errors = parse(raw)
    if RATES[1]*(times[-1]-times[0]) >= F(9, 5):
        return False
    for k in (0, len(times)-1):
        x, y = positions[k][0]-79, positions[k][1]+F(681, 10)
        error = errors[k]
        if x+error >= 0 or 6*(abs(y)+error) >= -(x+error):
            return False
    return True


def pairwise(intervals, times, rates=RATES):
    """Closed-form difference-constraint solution, independent of propagation."""
    if any(interval is None for interval in intervals):
        return None
    a, c = rates
    for j in range(len(times)):
        for i in range(j+1):
            dt = times[j]-times[i]
            if (intervals[i][0]+a*dt > intervals[j][1]
                    or intervals[j][0] > intervals[i][1]+c*dt):
                return None
    greatest = tuple(min(intervals[j][1]+(c if t >= times[j] else a)*(t-times[j])
                         for j in range(len(times))) for t in times)
    final = (max(interval[0]+a*(times[-1]-t) for interval, t in zip(intervals, times)),
             min(interval[1]+c*(times[-1]-t) for interval, t in zip(intervals, times)))
    return final, greatest


def decode(raw, method='DMI', steps=24):
    if method not in ('DC', 'MI', 'DMI') or type(steps) is not int or not 0 <= steps <= 64:
        raise ValueError('declared method and bounded integer inverse precision required')
    times, positions, errors = parse(raw)
    if not coverage_guard(raw):
        return dict(status='OUTSIDE_COVERAGE', selected=None, branches={}, evaluations=0, pruned=0)
    bounds = geometry_certificate()['derivative_bounds']
    duration = times[-1]-times[0]
    calls, pruned, branches = 0, 0, {}
    for name, domain in CHARTS.items():
        # Integrating the first-coordinate derivative enclosure gives a
        # necessary endpoint chord interval for every admitted timing law.
        directed = _interval_product(bounds[name][0], RATES)
        chord, error = positions[-1][0]-positions[0][0], errors[0]+errors[-1]
        if method == 'DMI' and not duration*directed[0]-error <= chord <= duration*directed[1]+error:
            branches[name] = dict(state='EMPTY', outer=None, witness=None)
            pruned += 1
            continue
        outers, inners = [], []
        for point, error in zip(positions, errors):
            outer = inner = domain
            for axis in (0, 1):
                if outer is None:
                    break
                sign = 1 if bounds[name][axis][0] > 0 else -1

                def evaluate(q):
                    nonlocal calls
                    calls += 1
                    return sign*curve(q)[axis]

                values = sorted((sign*(point[axis]-error), sign*(point[axis]+error)))
                out, inside = preimage(evaluate, values, domain, steps)
                outer, inner = intersect(outer, out), intersect(inner, inside)
            outers.append(outer)
            inners.append(inner)
        if method == 'DC':
            out, inside = pairwise(outers, times), pairwise(inners, times)
            final = out[0] if out is not None else None
            values = inside[1] if inside is not None else None
        else:
            out, inside = propagate(outers, times, RATES), propagate(inners, times, RATES)
            final = out[-1] if out is not None else None
            values = witness(inside, times, RATES) if inside is not None else None
        if values is not None and not verify_witness(values, raw, domain):
            raise ArithmeticError('constructed full-history witness failed verification')
        branches[name] = dict(state=('EMPTY' if final is None else
                                     'WITNESSED' if values is not None else 'POSSIBLE'),
                              outer=final, witness=values)
    retained = [name for name, branch in branches.items() if branch['state'] != 'EMPTY']
    if not retained:
        status, selected = 'INCOMPATIBLE', None
    elif len(retained) == 1 and branches[retained[0]]['state'] == 'WITNESSED':
        status, selected = 'UNIQUE', retained[0]
    elif all(branches[name]['state'] == 'WITNESSED' for name in retained):
        status, selected = 'AMBIGUOUS', None
    else:
        status, selected = 'UNRESOLVED', None
    return dict(status=status, selected=selected, branches=branches,
                evaluations=calls, pruned=pruned)
