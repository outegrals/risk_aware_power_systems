# Reviewer comments 2 and 4: risk interpretation and bidding limits

Reproduce the numerical evidence from the repository root:

```sh
python3 scripts/validate_risk_and_caps.py
```

The script uses Python's standard library and writes the sensitivity results, known capped cases, randomized checks, and `validation_report.md` in this directory.

- Comment 2: corrected the covariance expansion of mean-square price deviation; decomposed it into real-time price variance and squared price bias; derived covariance independence of bids with means, coefficients, and deterministic bounds fixed. The sensitivity sweep includes covariance offsets in expected profit. An independent 81-node Gaussian quadrature calculation checks the original settlement formula, moments, and unilateral objectives. The Gaussian distribution is not truncated or clipped. ERS increases all three operators' realized-profit standard deviations in this example; the paper does not claim variance or tail-risk reduction.
- Comment 4: proved a unique constrained equilibrium through a scalar residual equation; derived the constrained aggregate-profit optimum and a necessary-and-sufficient coefficient condition; constructed a feasible common-coefficient design. Validation separates equilibria at the old ERS coefficient from coordinated retuning. It includes the five-operator counterexample to the old midpoint partition, cap shortages, zero caps, and 200 equilibria across 100 heterogeneous instances.

The new Python script is the reference implementation for these revised tables. The older MATLAB scripts and saturated-branch figures are preserved as historical artifacts and are not evidence for the corrected constrained results. The manuscript no longer includes those old saturated plots.

These checks do not simulate physical microgrid dynamics, network congestion, communication failures, or cryptographic overhead. Forecast-relative caps do not guarantee realized delivery. Constrained controller design is outside the masked-consensus privacy guarantee.

Reviewed manuscripts: `main.tex` and `main_LONG.tex`. Compiled revision PDFs are saved in `output/pdf/` with suffix `_risk_caps_revision.pdf`. Original submitted and earlier revision PDFs are retained.
