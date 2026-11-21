# Native MATLAB revision-validation suite

Run `addpath('matlab'); run_validation('full')` from the repository root. Individual modules accept an output directory; dynamic and communication modules also accept `'full'` or `'smoke'`. Run dynamics before dispatch and plotting.

| Python predecessor | Native MATLAB entry point |
|---|---|
| scripts/validate_risk_and_caps.py | validate_risk_and_caps.m |
| scripts/validate_masked_consensus.py | validate_masked_consensus.m |
| scripts/validate_dynamic_microgrids.py | validate_dynamic_microgrids.m |
| scripts/audit_dynamic_microgrids.py | audit_dynamic_microgrids.m |
| scripts/validate_dispatch_feasibility.py | validate_dispatch_feasibility.m |
| scripts/validate_communication.py | validate_communication.m |
| scripts/render_microgrid_evidence.py | render_microgrid_evidence.m |
| output/incentives/check_contracted.py, skeptical_checks.py | validate_incentives.m |

Shared native numerical kernels are in `+rap`. There are no `py.*` calls, Python subprocesses, or Python execution dependencies.

## Reproduction versus new random experiments

`fixtures/microgrid_inputs.json` records original Python-generated asset sizes, phases, and discrete shock indices. It includes 30 replicates for each of n=10,50,100 and the two-operator audit assets. Shocks are reused across uncertainty scales as in the original paired study. MATLAB independently calculates every physical transition, conditional moment, cap, equilibrium, and cash flow. Every trajectory is compared with `output/microgrid/results.json` using relative tolerance 1e-7 (absolute tolerance 1e-7 for small values). The absolute maximum discrepancy is reported.

Privacy and communication replay public masks and input forecasts from the existing JSON evidence. These values are public test data, never operational cryptographic secrets. Communication iteration budgets are recomputed with upward-rounded Java BigDecimal arithmetic, not copied. Exact three-operator consensus uses Java BigInteger rational numerators/denominators. The executable exhaustive p=7 privacy check compares full multiset views.

Random capped cases, incentive cases, and packet-loss episodes use MATLAB's seeded Twister generator. They are fresh checks under the same model, not bit-for-bit NumPy/Python random streams. Deterministic analytical values, inequalities and target recovery are checked directly. Do not replace original reported random-sample numbers with new ones without reviewing provenance.

## Outputs and limitations

- Independent microgrid audit enumerates 81 joint outcomes at each of 80 two-operator states, evaluates settlements directly, and fits unilateral objectives. It compares sorted roots against an independent bisection solver.
- Dispatch checks cover the full current-interval support at visited physical states and the entire bid box, not all possible multi-interval histories.
- Gaussian 81-node quadrature checks expected revenue, revenue variance, squared price spreads, and directly fitted best responses. A separate bounded 16-state distribution preserves first/second moments, not Gaussian fourth moments.
- Communication runs all 32 original topology/size/precision profiles. Counts above 12,000 use spectral evaluation rather than pretending every iteration was executed. Kernel timing measures at most 2,000 actual updates per trial; the count is explicitly reported. Timings are MATLAB/host-specific and do not reproduce Python host timings.
- Buffered retries preserve synchronous lossless states; synthetic network times assume reliable retries and independent attempt loss. One-peer scheduling uses graph edge coloring. Max-geometric tails are bounded; multinomial sampling uses inverse-binomial CDFs via base MATLAB `betainc`.
- Probability intervals are unadjusted Wilson intervals; dynamic cash intervals use independent trajectory replicates. Smoke-mode dynamic confidence intervals are omitted/null because one replicate is insufficient.
- Incentive analyses remain exploratory and are not incorporated into the original working manuscript. Their fixed-input, binding-settlement, and participation limitations remain in `output/incentives`.
- MATLAB exports native vector PDF/PNG trajectory graphics, a TikZ fragment, and a convergence plot. These are new artifacts; old manuscript graphics are not overwritten.
- JSON schemas are organized for MATLAB and are not byte-for-byte copies of historical schemas. Original evidence is retained.

## Why legacy scripts are retained

Root power-market scripts remain unchanged for historical reproducibility and manuscript figure paths. The new runner is the central validation workflow. Unrelated MATLAB exercises are archived separately, with a recovery manifest.
