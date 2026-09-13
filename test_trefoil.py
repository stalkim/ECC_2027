"""Focused checks for the sourced geometry and exact local decoder."""
import copy
from fractions import Fraction as F
from itertools import product
from math import atan, cos, sin
import unittest

from straight_history import propagate
from trefoil_benchmark import endpoints, prepare, prepare_case
from trefoil_history import (CHARTS, LOCAL_PRIOR, RATES, coverage_guard, curve, decode,
                             expanded_curve, geometry_certificate, pairwise,
                             verify_geometry_certificate, verify_witness)


class TrefoilTests(unittest.TestCase):
    def test_source_geometry_and_expanded_identity(self):
        for q in (F(-11, 10), F(-1), F(-9, 10), F(0), F(9, 10), F(1), F(11, 10)):
            self.assertEqual(curve(q), expanded_curve(q))
            u = 2*atan(float(q))
            radius = 160+80*cos(3*u)
            reference = (radius*cos(2*u)+79, radius*sin(2*u)-68.1)
            for exact, approximate in zip(curve(q), reference):
                self.assertAlmostEqual(float(exact), approximate, places=10)
        self.assertEqual(curve(1), curve(-1))
        self.assertEqual(curve(1), (F(-81), F(-681, 10)))

    def test_rigorous_geometry_certificate(self):
        certificate = geometry_certificate()
        self.assertTrue(verify_geometry_certificate(certificate))
        self.assertGreater(certificate['derivative_bounds']['negative'][0][0], 0)
        self.assertLess(certificate['derivative_bounds']['positive'][0][1], 0)
        for name in CHARTS:
            self.assertLess(certificate['derivative_bounds'][name][1][1], 0)
        corrupted = copy.deepcopy(certificate)
        corrupted['derivative_bounds']['negative'] = ((F(250), F(251)), (F(-366), F(-237)))
        self.assertFalse(verify_geometry_certificate(corrupted))

    def test_two_chart_diagnostic(self):
        for name in CHARTS:
            case, truth = prepare_case(name, 4, F(4, 5), 'zero')
            raw = case['history']
            self.assertTrue(verify_witness(truth['parameters'], raw, CHARTS[name]))
            results = [decode(raw, method) for method in ('DC', 'MI', 'DMI')]
            self.assertEqual(len({r['status'] for r in results}), 1)
            self.assertNotIn(results[0]['status'], ('INCOMPATIBLE', 'UNRESOLVED', 'OUTSIDE_COVERAGE'))
            for result in results:
                self.assertEqual(set(result['branches']), set(CHARTS))
                outer = result['branches'][name]['outer']
                self.assertLessEqual(outer[0], truth['parameters'][-1])
                self.assertGreaterEqual(outer[1], truth['parameters'][-1])
                for branch, data in result['branches'].items():
                    if data['witness'] is not None:
                        self.assertTrue(verify_witness(data['witness'], raw, CHARTS[branch]))
            for branch in CHARTS:
                self.assertEqual(results[0]['branches'][branch]['outer'], results[1]['branches'][branch]['outer'])

    def test_boundary_status_is_not_guessed(self):
        raw = dict(prior=LOCAL_PRIOR, times=(F(0),), positions=(curve(F(1)),), linf_errors=(F(0),))
        self.assertEqual(decode(raw)['status'], 'AMBIGUOUS')
        raw['positions'] = (curve(F(101, 100)),)
        self.assertEqual(decode(raw, steps=0)['status'], 'UNRESOLVED')
        # An impossible interior record retains covered endpoints.
        raw['times'] = (F(0), F(1, 100), F(1, 50))
        raw['positions'] = (curve(1), (F(10000), F(10000)), curve(1))
        raw['linf_errors'] = (F(1), F(1), F(1))
        self.assertEqual(decode(raw)['status'], 'INCOMPATIBLE')

    def test_explicit_local_prior_and_invalid_inputs(self):
        case, _ = prepare_case('negative', 4, F(4, 5), 'zero')
        raw = copy.deepcopy(case['history'])
        del raw['prior']
        self.assertEqual(decode(raw)['status'], 'AMBIGUOUS')
        raw['prior'] = LOCAL_PRIOR
        for key, value in (('times', (0, 0, 1, 2, 3)), ('positions', [(0, 0, 0)]*5),
                           ('linf_errors', [-1]*5), ('times', ())):
            invalid = dict(raw, **{key: value})
            with self.assertRaises(ValueError):
                decode(invalid)
        for method, steps in (('other', 24), ('MI', -1), ('MI', 65), ('MI', True)):
            with self.assertRaises(ValueError):
                decode(raw, method, steps)

    def test_coverage_guard_cannot_be_bypassed_and_rejects_equality(self):
        raw = dict(prior=LOCAL_PRIOR, times=(F(0),), positions=((F(0), F(0)),),
                   linf_errors=(F(0),))
        self.assertFalse(coverage_guard(raw))
        self.assertEqual(decode(raw)['status'], 'OUTSIDE_COVERAGE')
        # X=-6,Y=1 is exactly on the cone boundary; strict coverage rejects it.
        raw['positions'] = ((F(73), F(1)-F(681, 10)),)
        self.assertFalse(coverage_guard(raw))
        self.assertEqual(decode(raw)['status'], 'OUTSIDE_COVERAGE')
        raw['positions'] = (curve(1), curve(1))
        raw['linf_errors'] = (F(0), F(0))
        raw['times'] = (F(0), F(9, 5)/RATES[1])
        self.assertFalse(coverage_guard(raw))
        self.assertEqual(decode(raw)['status'], 'OUTSIDE_COVERAGE')

    def test_pairwise_against_propagation(self):
        times = (F(0), F(1), F(3))
        candidates = (None, (F(0), F(1, 20)), (F(1, 20), F(1, 10)), (F(0), F(1, 5)))
        for intervals in product(candidates, repeat=3):
            direct, recursive = pairwise(intervals, times), propagate(intervals, times, RATES)
            self.assertEqual(direct is None, recursive is None)
            if direct is not None:
                self.assertEqual(direct[0], recursive[-1])
                q = direct[1]
                self.assertTrue(all(a <= v <= b for (a, b), v in zip(intervals, q)))
                self.assertTrue(all(RATES[0]*(v-u) <= b-a <= RATES[1]*(v-u)
                                    for u, v, a, b in zip(times, times[1:], q, q[1:])))

    def test_fixed_protocol_and_truth_separation(self):
        cases, truth = prepare()
        self.assertEqual(len(cases), 36)
        self.assertEqual(len(truth), 36)
        for case in cases:
            raw, actual = case['history'], truth[case['case']]
            self.assertEqual(set(raw), {'prior', 'times', 'positions', 'linf_errors'})
            self.assertTrue(verify_witness(actual['parameters'], raw, CHARTS[actual['branch']]))
            self.assertTrue(coverage_guard(raw))
            self.assertEqual(raw['times'][-1]-raw['times'][0], F(case['intervals'], 50))
            self.assertEqual(len(endpoints(raw)['times']), 2)
        case, actual = prepare_case('positive', 4, F(4, 5), 'zero')
        invalid = tuple(q+1 for q in actual['parameters'])
        self.assertFalse(verify_witness(invalid, case['history'], CHARTS['positive']))


if __name__ == '__main__':
    unittest.main()
