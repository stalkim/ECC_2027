# ECC 2027: Certified Branch Recovery

Code for **Sharp Distinguishability Bounds and Certified Branch Recovery on
Self-Intersecting Paths**, by Stanislav Kim and Anton Pyrkin.

The decoder finds which branches of a known path can explain a finite
position history under bounded errors and an unknown positive progress rate.

## Run

Python 3.10 or later; no third-party packages. From the repository root:

```bash
python3 demo.py
python3 -m unittest discover -s . -p 'test_*.py'
python3 trefoil_benchmark.py --output results/trefoil-01
```

On Windows, use `python` if needed. Run without `-O`, since the retained
regression checks use assertions. Each comparison requires a new output
directory and refuses to overwrite an existing one.

The demo returns ambiguous endpoints and a unique full history, with no
directional pruning. It uses one case from the comparison below.

## Source and assumptions

The path is the trefoil from Yao et al., *Singularity-Free Guiding Vector
Field for Robot Navigation*, IEEE T-RO 37(4), 1206–1221 (2021),
[Section VI-D](https://arxiv.org/pdf/2012.01826v3),
[doi:10.1109/TRO.2020.3043690](https://doi.org/10.1109/TRO.2020.3043690).
Its sourced geometry, in metres, is
`(80 cos(3u)+160)(cos(2u),sin(2u))+(79,-68.10)`.

The exact substitution `q=tan(u/2)` permits rational coordinate evaluations.
The two crossing charts are `[-1.1,-0.9]` and `[0.9,1.1]`. A position-based
endpoint guard verifies their coverage; all methods consider both branches.
Coordinate monotonicity is checked with rational interval bounds, not samples.

The comparison is a kinematic adaptation, not a replay of flight logs or the
source controller. Source guidance frequency motivates 50 Hz position
sampling. The 12 m/s reference and 400 m-per-q crossing tangent give nominal
`q_dot=0.03`; `[0.027,0.033]` is our imposed 10% rate band, not a measured
ground-speed bound. The error bounds `0.8,1.6,4` m are 1%, 2%, and 5% of
the sourced 80 m scale; they are sensitivity settings, not sensor calibration.

## Comparison

All 36 cases combine two crossing passages, spans `0.08,0.32,1.28` s,
three error bounds, and zero/alternating bounded errors. Alternating errors
use 80% of the allowed first-coordinate error. Cases share six underlying
motions. There is no claim of 36 independent flights or a statistical test.

| Method | UNIQUE | AMBIGUOUS | Coordinate calls |
| --- | ---: | ---: | ---: |
| Endpoint MI | 20 | 16 | 15552 |
| Full DC | 26 | 10 | 225504 |
| Full MI | 26 | 10 | 225504 |
| Full DMI | 26 | 10 | 126576 |

Six cases gain uniqueness from intermediate samples without pruning.
DMI saves 43.9% of coordinate calls over full MI; this is not a runtime claim.
DC checks pairwise temporal constraints; MI propagates feasible intervals;
DMI adds lossless directional pruning. All share the inverse builder with
24 bisections. Witnesses are checked in an expanded rational path equation.
DC/MI current-parameter outers agree. DC is not a COMMA implementation.

The output directory contains exact inputs, separate truth, all decisions,
witnesses, source hashes, derivative certificates, summaries and CSV files.
Truth is never passed to the decoder. `UNRESOLVED` is preserved if finite
precision cannot establish feasibility or exclusion. `INCOMPATIBLE` means
no branch fits; `OUTSIDE_COVERAGE` means the local guard is not established.
Neither outcome diagnoses a sensor fault.

`trefoil_history.py` and `trefoil_benchmark.py` implement this comparison.
`monotone_history.py` and `straight_history.py` provide shared inverse,
interval-propagation and witness routines. The two test modules check these
routines and the trefoil implementation. Earlier experiments remain in Git
history. The repository excludes the manuscript, flight dynamics, global
lap recovery and closed-loop control.
