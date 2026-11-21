# Risk-aware power-system bidding research

## Manuscripts

- `main_before_reviewer_changes.tex`: recovered original, current working manuscript for targeted reviewer changes.
- `main.tex`: earlier revision branch of the manuscript, retained for reference.
- `main_LONG.tex`: expanded manuscript; do not modify without a new explicit request.
- `output/reviewer_notes/all_reviewer_issues.md`: consolidated issue register.

Code experiments do not imply approval to incorporate a new model into the manuscript.

## MATLAB validation (current workflow)

MATLAB R2025b, base MATLAB plus its bundled Java runtime. No Python runtime or additional MATLAB toolboxes are required.

From the repository root:

```matlab
addpath('matlab');
run_validation('full');  % all 360 trajectories and all communication profiles
% run_validation('smoke'); % one trajectory per setting, small communication graphs
```

Results go to `output/matlab/full` or `output/matlab/smoke`. Original published Python results are preserved in their existing output folders. The runner never edits a manuscript, bibliography, or original evidence artifact.

See [matlab/README.md](matlab/README.md) for the source mapping, input provenance, scope, and parity checks.

## Existing MATLAB work

The original power-market scripts remain at the root to preserve their relative paths and historical behavior:
`distributed_risk_aware_nash.m`, `risk_aware_power_system.m`, `stochastic_risk_aware_nash_mc.m`, `central_pwr_sys.m`, `distributed_pwr_sys.m`, `consensus_test.m`, and `test.m`.

These are historical model implementations, not substitutes for the revised validation suite. In particular, do not infer a privacy guarantee or size-independent consensus rate from their comments. `test.m` is retained because it is a power-market prototype, not an unrelated test.

## Archived material

Unrelated UAV/encirclement, formation homework, switching-system exercises, and banner code have been moved intact to `archive/unrelated_matlab`. `manifest.json` records original paths and SHA-256 hashes. Nothing was permanently deleted.

Python predecessors are archived under `archive/python_reference` after MATLAB verification. Their previously generated artifacts retain their original provenance. Archived scripts are historical source snapshots; their old relative paths need restoration before rerunning them. They are not called by MATLAB.
