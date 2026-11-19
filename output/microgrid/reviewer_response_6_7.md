# Draft response to Reviewer 1, comments 6 and 7

This is a local draft for the response letter. It has not been submitted.

## Comment 6: Dynamic and cyber-physical microgrid connection

Thank you for pointing out the gap between the original market formulation and its microgrid motivation. We have revised the title to “A Distributed Risk-Aware Bidding Strategy for Renewable Operators and Microgrid Aggregators” and explicitly placed the contribution at the supervisory power-exchange layer. The revised formulation adds battery energy-state dynamics, charging/discharging efficiency and power bounds, local renewable generation and demand, feeder export limits, curtailment, and grid imports. It distinguishes that energy-management application from voltage/frequency control and electrical network feasibility.

We added a 48-hour sequential microgrid study with 15-minute market intervals. A prescribed local battery policy operates independently of market bids. Its current energy state determines the conditional export distribution and an availability-based bid cap before the next disturbance is sampled. This independence is necessary for the existing interval-level bidding analysis: imports and battery-throughput costs contribute bid-independent terms, while the bid-dependent profit expression and capped best responses remain valid. We do not claim that the theorem solves optimal storage scheduling, establishes a multi-period strategic equilibrium, or covers a battery policy that follows the chosen bid. The sequential study updates one-interval commitments; it is not a fixed 24-hour schedule submitted on the preceding day.

We explicitly report that all bidding strategies share their physical trajectories and terminal battery energy. Differences in scheduled-export shortfalls concern market commitments, not improved voltage/frequency regulation or reduced unserved local demand. The study uses synthetic, export-dominated profiles, unrestricted import support, and no terminal-energy valuation. These limits are stated in the revised paper and reproducibility materials.

## Comment 7: Numerical validation

We substantially expanded the numerical evidence in both manuscripts. The earlier revision added a fixed-mean/fixed-cap uncertainty sweep, independent Gaussian quadrature checks, and 200 capped equilibria across 100 heterogeneous instances. The present revision additionally includes:

- 360 synthetic 48-hour microgrid trajectories: 10, 50, and 100 operators; four disturbance amplitudes; and 30 independent replicates per setting. Each replicate pairs RN, fixed ERS, and coordinated retuning on identical physical trajectories. This covers 69,120 market intervals and 207,360 bidding equilibria. The paper reports net cash including import purchases and throughput cost, scheduled-export shortfalls, and approximate paired confidence intervals across replicates. Results do not imply a universal improvement in all risk measures: retuning slightly increases shortfalls relative to fixed ERS in the reported examples.
- Physical and equilibrium checks across the declared conditional support. There were no state-of-charge or power-limit violations; maximum power-balance and normalized best-response residuals were below 4.5e-16 MW and 3.9e-16, respectively. A separate implementation enumerates 6,480 joint scenarios, checking settlement expectations, directly fitted unilateral best responses, availability probabilities, and an independent bisection solver.
- A separate uncapped masked-consensus study with 3, 10, 50, 100, and 200 operators; star, ring, degree-at-most-four circulant, and complete graphs; two encoding grids; and 96 public-mask trials across 32 cases. It reports conservative spectral and binary64 decoding bounds, logical message counts, observed errors, and one-host computation measurements. Long-ring numerical diagnostics evaluate ideal iterates spectrally; timing runs execute at most 2,000 updates, so we do not present those timings as full distributed runs.
- Synthetic network timing under 1/10 ms attempt slots, independent failure probabilities 0/0.05/0.20, buffered reliable retries, and either parallel links or one-peer-per-node scheduling. An event-level replay checks that retries preserve the lossless logical states. The model charges mask exchanges and exposes hub degree. At 100 operators, 10 ms attempts, and 20% failure with one-peer scheduling, the star and ring miss the 15-minute deadline in every one of 2,000 sampled episodes, whereas the degree-four and complete graphs meet it in every episode. These are modeled outcomes, not deployed-network guarantees.

The communication study concerns uncapped Algorithm 1; it does not supply a private constrained implementation for the physical reference. No physical fallback is simulated after a missed communication deadline, and cryptographic processing, synchronization/barrier traffic, correlated losses, and participant dropout remain outside the timing model. The short paper stays within the 10-page limit; the expanded paper contains full profiles, accounting, arithmetic bounds, and additional interpretation.

## Artifacts

- `scripts/validate_dynamic_microgrids.py`, `scripts/audit_dynamic_microgrids.py`, and `scripts/render_microgrid_evidence.py`.
- `output/microgrid/results.json`, `independent_audit.json`, `representative_trajectory.csv`, `README.md`, and `validation_report.md`.
- `scripts/validate_communication.py` and `output/communication/communication_validation.json`, `README.md`, and `convergence_n100.svg`.

These changes provide a substantially stronger, reproducible supervisory application and communication study. They do not demonstrate a complete cyber-physical deployment or establish that the reviewer will consider either comment closed.
