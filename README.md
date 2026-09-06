# ACC 2027: Certified Branch Recovery

Code for **Sharp Distinguishability Bounds and Certified Branch Recovery on
Self-Intersecting Paths**, by Stanislav Kim and Anton Pyrkin.

The code determines which directed branches can explain a finite position
history when the path is known, the speed profile is unknown, and measurement
errors and progress rates have known bounds.

## Requirements

Python 3.10 or later. No third-party packages are needed. Run commands from
the repository root without `-O`, because validation uses assertions.

```bash
python3 demo.py
python3 -m unittest discover -s . -p 'test_*.py'
```

On Windows, use `python` instead of `python3` if needed.

## Reproduce the results

Each command requires a new output directory and will not overwrite an
existing one.

Geometric boundary, curvature, informative histories, and scaling:

```bash
python3 geometric_checks.py --output results/geometric-01
```

Expected results:

- 48 ambiguous midpoint histories at or below the straight-branch threshold
  and 32 incompatible midpoint histories above it;
- all 256 model-valid histories above the threshold uniquely identify the
  correct branch;
- nine cases in which endpoints are ambiguous but the full history is unique;
- 369 chart-covered circular-arc cases, with 296 certified by the general
  bound and 305 by the projected-chord bound;
- unique decisions without pruning for 3, 9, and 33 samples.

Polynomial-path comparison and precision check:

```bash
python3 benchmark.py --check-only --output results/polynomial-check-01
python3 check_precision.py results/polynomial-check-01
python3 benchmark.py --output results/polynomial-timing-01
```

The 60 fixed histories contain 48 `UNIQUE`, 6 `AMBIGUOUS`, and 6
`INCOMPATIBLE` cases. DC and MI use 57024 coordinate evaluations per batch;
DMI uses 28512. DC checks pairwise temporal constraints, MI propagates
feasible intervals, and DMI adds directional pruning. All three methods share
the coordinate inverse routine.

The scripts save inputs, separate ground truth, decisions, witnesses, source
hashes, timings, and summaries under `results/`. Ground truth is used only for
validation and is not passed to the decoders. The checks are deterministic.

Reference timings in the paper used Python 3.10.12 on Linux x86-64, one Python
process, and an Intel Xeon Platinum 8468 CPU. Batch timings exclude setup and
file output and will vary by machine.

## Decisions

- `UNIQUE`: exactly one branch has a checked feasible history.
- `AMBIGUOUS`: both branches have checked feasible histories.
- `INCOMPATIBLE`: neither branch can explain the measurements.
- `UNRESOLVED`: finite-precision bounds do not settle the decision.
- `OUTSIDE_COVERAGE`: the required chart-coverage condition fails.

A nonempty outer interval alone does not prove that a feasible history exists.

## Files

- `straight_history.py`: straight branches with Euclidean error balls;
- `monotone_history.py`: interval propagation and cubic-path decoder;
- `extension_history.py`: quartic/quintic paths and DC, MI, DMI;
- `geometric_checks.py`: geometric and informative-history checks;
- `benchmark.py` and `check_precision.py`: polynomial comparison;
- `demo.py` and `test_*.py`: examples and regression tests.

## Scope

The repository covers the stated known-path models. It does not include the
manuscript, arbitrary-curve recovery, unknown geometry, global lap recovery,
vehicle dynamics, or closed-loop control.

The temporal-constraint formulation follows R. Dechter, I. Meiri, and
J. Pearl, ["Temporal constraint
networks"](https://www.sciencedirect.com/science/article/pii/0004370291900066),
*Artificial Intelligence*, 49 (1991), 61-95. Interval reachability and
set-membership estimation are established methods; DC is not an external
map-matching implementation.
