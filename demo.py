"""A short local recovery example on the published trefoil geometry."""
from fractions import Fraction as F

from trefoil_benchmark import endpoints, prepare_case
from trefoil_history import decode


def main():
    # First duration/error setting in the fixed sensitivity grid; no tuning.
    case, _ = prepare_case('positive', 4, F(4, 5), 'alternating')
    print('Published trefoil geometry; synthetic bounded-error positions.')
    print('50 Hz, span 0.08 s, coordinate-error bound 0.8 m.')
    for name, history in (('endpoints', endpoints(case['history'])),
                          ('full history', case['history'])):
        result = decode(history, 'DMI')
        branch = result['selected'] or 'none'
        print(f'{name}: {result["status"]}; branch: {branch}; '
              f'pruned: {result["pruned"]}')


if __name__ == '__main__':
    main()
