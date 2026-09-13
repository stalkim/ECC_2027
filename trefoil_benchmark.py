"""Run the fixed 36-case, source-derived local trefoil comparison once.

The source geometry and 50 Hz update frequency come from Yao et al. (2021).
The constant parameter rate, rate band and error sensitivities are explicitly
declared kinematic adaptations, not a replay of source flight measurements.
The decoder checks chart coverage from the measured endpoint boxes and rate
cap; the optional legacy prior field is metadata and provides no admission.
"""
import argparse
from collections import Counter
import csv
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import platform
import sys

from trefoil_history import (CHARTS, LOCAL_PRIOR, RATES, SOURCE, decode,
                             expanded_curve, geometry_certificate,
                             verify_geometry_certificate, verify_witness)

PRECISION = 24
REFERENCE_RATE = F(3, 100)


def prepare_case(branch, intervals, error, profile):
    if branch not in CHARTS or intervals not in (4, 16, 64):
        raise ValueError('undeclared chart or sample schedule')
    error = F(error)
    if error not in (F(4, 5), F(8, 5), F(4)) or profile not in ('zero', 'alternating'):
        raise ValueError('undeclared sensitivity setting')
    center = F(-1 if branch == 'negative' else 1)
    times = tuple((F(m)-F(intervals, 2))/50 for m in range(intervals+1))
    parameters = tuple(center+REFERENCE_RATE*t for t in times)
    positions = []
    for m, q in enumerate(parameters):
        x, y = expanded_curve(q)
        if profile == 'alternating':
            x += F(4, 5)*error*(-1)**m
        positions.append((x, y))
    identifier = f'{branch}_N{intervals}_b{str(error).replace("/", "_")}_{profile}'
    raw = dict(prior=LOCAL_PRIOR, times=times, positions=positions,
               linf_errors=(error,)*len(times))
    case = dict(case=identifier, branch_center=str(center), intervals=intervals,
                error=error, profile=profile, history=raw)
    truth = dict(branch=branch, parameters=parameters)
    return case, truth


def prepare(*, diagnostic=False):
    settings = ((branch, 4, F(4, 5), 'zero') for branch in CHARTS) if diagnostic else product(
        CHARTS, (4, 16, 64), (F(4, 5), F(8, 5), F(4)), ('zero', 'alternating'))
    cases, truth = [], {}
    for settings_row in settings:
        case, actual = prepare_case(*settings_row)
        cases.append(case)
        truth[case['case']] = actual
    return cases, truth


def endpoints(raw):
    selected = {key: (raw[key][0], raw[key][-1])
                for key in ('times', 'positions', 'linf_errors')}
    if 'prior' in raw:
        selected['prior'] = raw['prior']
    return selected


def _check(result, raw, truth, *, endpoint_only=False):
    if result['status'] in ('INCOMPATIBLE', 'OUTSIDE_COVERAGE'):
        raise ArithmeticError('a valid history lost all model-compatible branches')
    if set(result['branches']) != set(CHARTS):
        raise ArithmeticError('decoder did not account for every supplied chart')
    if result['selected'] is not None and result['selected'] != truth['branch']:
        raise ArithmeticError('incorrect unique branch')
    actual = (truth['parameters'][0], truth['parameters'][-1]) if endpoint_only else truth['parameters']
    if not verify_witness(actual, raw, CHARTS[truth['branch']]):
        raise ArithmeticError('generator violated geometry, chart, error or rate assumptions')
    outer = result['branches'][truth['branch']]['outer']
    if outer is None or not outer[0] <= actual[-1] <= outer[1]:
        raise ArithmeticError('true current parameter absent from returned outer interval')
    witnesses = 0
    for name, branch in result['branches'].items():
        if branch['witness'] is not None:
            if not verify_witness(branch['witness'], raw, CHARTS[name]):
                raise ArithmeticError('full-history witness failed direct equation check')
            if branch['outer'] is None or not branch['outer'][0] <= branch['witness'][-1] <= branch['outer'][1]:
                raise ArithmeticError('witness current parameter absent from outer interval')
            witnesses += 1
    return witnesses


def _json(path, value):
    path.write_text(json.dumps(value, default=str, indent=2, sort_keys=True)+'\n')


def _digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run(output, *, diagnostic=False):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    directory = Path(__file__).resolve().parent
    paths = [directory/name for name in ('trefoil_history.py', 'trefoil_benchmark.py',
                                         'test_trefoil.py', 'monotone_history.py',
                                         'straight_history.py')]
    protocol = directory.parents[1]/'research/acc2027_sourced_validation_20260913.md'
    if protocol.exists():
        paths.append(protocol)
    pins = {p.name: _digest(p) for p in paths}
    cases, truth = prepare(diagnostic=diagnostic)
    _json(output/'inputs.json', cases)
    _json(output/'truth.json', truth)
    certificate = geometry_certificate()
    _json(output/'geometry_certificate.json', certificate)
    _json(output/'manifest.json', dict(
        source=SOURCE, source_section='VI-D, trefoil equation, physical coordinates in metres',
        source_geometry='(80 cos(3u)+160)(cos(2u),sin(2u))+(79,-68.10)',
        source_parameter_relation='u=0.0045*w_source; q=tan(u/2)',
        source_frequency_hz=50, local_charts=CHARTS, local_prior=LOCAL_PRIOR,
        chart_coverage='endpoint error-box cone X<0, |Y/(-X)|<1/6 and c*T<1.8',
        reference_parameter_rate=REFERENCE_RATE, allowed_parameter_rates=RATES,
        error_settings_metres=(F(4, 5), F(8, 5), F(4)),
        adaptation='synthetic constant-parameter-rate histories; no flight-data/controller replay',
        inverse_bisections=PRECISION, diagnostic=diagnostic, cases=len(cases),
        source_hashes=pins, python=sys.version, platform=platform.platform(),
        inputs_sha256=_digest(output/'inputs.json'), truth_sha256=_digest(output/'truth.json')))
    rows, comparisons = [], []
    try:
        if not verify_geometry_certificate(certificate):
            raise ArithmeticError('geometry certificate failed rational verification')
        for case in cases:
            raw, actual = case['history'], truth[case['case']]
            if not verify_witness(actual['parameters'], raw, CHARTS[actual['branch']]):
                raise ArithmeticError('source-derived nominal history failed admission')
            outputs = {}
            for label, method, data in (('endpoint_MI', 'MI', endpoints(raw)),
                                        ('DC', 'DC', raw), ('MI', 'MI', raw), ('DMI', 'DMI', raw)):
                result = decode(data, method, steps=PRECISION)
                checked = _check(result, data, actual, endpoint_only=(label == 'endpoint_MI'))
                row = dict(case=case['case'], method=label, result=result,
                           checked_witnesses=checked, truth_in_current_outer=True)
                rows.append(row)
                with (output/'calls.jsonl').open('a') as handle:
                    handle.write(json.dumps(row, default=str, sort_keys=True)+'\n')
                outputs[label] = result
            if any(outputs['DC']['branches'][n]['outer'] != outputs['MI']['branches'][n]['outer']
                   for n in CHARTS):
                raise ArithmeticError('DC/MI current-outer disagreement')
            if len({(outputs[m]['status'], outputs[m]['selected']) for m in ('DC', 'MI', 'DMI')}) != 1:
                raise ArithmeticError('full-history method status disagreement')
            entry = dict(case=case['case'], branch=actual['branch'], intervals=case['intervals'],
                         duration=float(raw['times'][-1]-raw['times'][0]),
                         error=float(case['error']), profile=case['profile'])
            for label, result in outputs.items():
                entry[label+'_status'] = result['status']
                entry[label+'_evaluations'] = result['evaluations']
                entry[label+'_pruned'] = result['pruned']
            comparisons.append(entry)
            if (actual['branch'], case['intervals'], case['error'], case['profile']) == (
                    'positive', 16, F(8, 5), 'alternating'):
                # This representative is fixed before decoding: the central
                # duration/error setting on the positive chart.
                _json(output/'representative.json', dict(case=case, truth=actual, results=outputs))
                with (output/'representative.csv').open('w', newline='') as handle:
                    writer = csv.writer(handle)
                    writer.writerow(('time_s', 'true_parameter', 'true_x_m', 'true_y_m',
                                     'observed_x_m', 'observed_y_m', 'coordinate_error_m'))
                    for t, q, point, error in zip(raw['times'], actual['parameters'],
                                                  raw['positions'], raw['linf_errors']):
                        x, y = expanded_curve(q)
                        writer.writerow(tuple(float(v) for v in (t, q, x, y, *point, error)))
            print(f'{case["case"]}: {outputs["endpoint_MI"]["status"]} -> '
                  f'{outputs["MI"]["status"]}; DMI pruned={outputs["DMI"]["pruned"]}', flush=True)
        if pins != {p.name: _digest(p) for p in paths}:
            raise ArithmeticError('source files changed during the run')
    except Exception as error:
        _json(output/'summary.json', dict(status='FAILED', cases_completed=len(comparisons),
                                         calls_completed=len(rows), error=str(error)))
        raise
    methods = {}
    for label in ('endpoint_MI', 'DC', 'MI', 'DMI'):
        chosen = [row for row in rows if row['method'] == label]
        methods[label] = dict(statuses=dict(Counter(row['result']['status'] for row in chosen)),
                              coordinate_evaluations=sum(row['result']['evaluations'] for row in chosen),
                              pruned_branches=sum(row['result']['pruned'] for row in chosen),
                              checked_witnesses=sum(row['checked_witnesses'] for row in chosen))
    transitions = Counter(row['endpoint_MI_status']+' -> '+row['MI_status'] for row in comparisons)
    summary = dict(status='COMPLETE', diagnostic=diagnostic, cases=len(cases), calls=len(rows),
                   all_truth_inclusions_verified=True, all_witnesses_checked=True,
                   dc_mi_current_outers_equal=True, full_statuses_equal=True,
                   source_hashes_unchanged=True, endpoint_to_full=dict(transitions), methods=methods)
    _json(output/'summary.json', summary)
    with (output/'case_summary.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparisons[0]))
        writer.writeheader()
        writer.writerows(comparisons)
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='new output directory; never overwritten')
    parser.add_argument('--diagnostic', action='store_true', help='only the two prescribed on-path diagnostics')
    arguments = parser.parse_args()
    run(arguments.output.resolve(), diagnostic=arguments.diagnostic)
