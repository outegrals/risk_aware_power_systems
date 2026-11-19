# Communication and convergence validation

Run `python3 scripts/validate_communication.py`. Python 3 and NumPy are sufficient.

These are mathematical and numerical checks plus a synthetic network model. The timing measurements are from one host running all nodes. They are not deployed microgrid latency, secure-channel measurements, or evidence of physical closed-loop stability.

## Protocol and guarantees

- The public consensus graph contains n operators and the utility. The private mask graph is an operator cycle; its honest-induced connectivity still limits collusion tolerance. Its private logical links are assumed directly available with the same attempt-slot model; additional relay/routing costs are excluded.
- The fine-grid profile uses public B=10000 integer units and Delta=.001 MW (forecasts at most 10 MW), at n=10,50,100. Its forecasts are independently generated; this is not a communication replay of the separate physical trajectory. A baseline profile uses B=120, Delta=1 MW at n=3,10,50,100,200. Each uses p>nB and fresh seeded public test masks. The script records inputs, masks, modular wraps, and decoded sums. These masks are not cryptographic randomness.
- W=I-L/(dmax+1). Star and complete spectra are exact; ring/circulant numerical spectra are accompanied by conservative rational bounds derived from sin(pi*t)>=2t for 0<=t<=1/2.
- The fixed public budget uses m*ceil(sqrt(m))*(p-1)*rho_upper^k<1/4. Upward-rounded decimal arithmetic verifies the bound. A conservative IEEE binary64 error bound also keeps the total below 1/2 for every case.
- The arithmetic certificate assumes correctly rounded binary64 addition, multiplication and division, round to nearest, exact small integer inputs/coefficients, and no overflow/underflow. It is not a hardware or compiler formal verification.
- First-correct-decoding rounds are offline oracle diagnostics. They are never used for input-dependent protocol stopping.

- Budgets above 12000 updates are evaluated through the known spectral form of W^k, not by executing every update. Other convergence cases execute all updates. Timing runs execute min(k,2000) updates, with their count reported explicitly; no extrapolated full-protocol runtime is claimed. Spectral evaluation is a fast numerical diagnostic, while the separate iteration/arithmetic bounds apply to the actual positive-sum kernel.

## Fixed-budget convergence

The table gives the maximum final error across three fresh public-mask trials. Directed messages exclude mask exchange, authentication, acknowledgments, and barrier control.

| B; Delta MW | Operators | Graph | Max degree | Numerical rho | Certified rho upper | Public rounds | Offline first correct, range | Final scaled error | Directed messages | Timed updates | Median host kernel ms |
|:---|---:|:---|---:|---:|---:|---:|:---|---:|---:|---:|---:|
| 120; 1.0 | 3 | star | 3 | 0.7500000 | 0.7500000 | 37 | 26–30 | 0.051458 | 222 | 37 | 0.076 |
| 120; 1.0 | 3 | ring | 2 | 0.3333333 | 0.6666667 | 27 | 8–8 | 2.7831e-10 | 216 | 27 | 0.131 |
| 120; 1.0 | 3 | circulant4 | 3 | 0.0000000 | 0.5000000 | 16 | 1–1 | 0 | 192 | 16 | 0.112 |
| 120; 1.0 | 3 | complete | 3 | 0.0000000 | 0.0000000 | 1 | 1–1 | 0 | 12 | 1 | 0.004 |
| 120; 1.0 | 10 | star | 10 | 0.9090909 | 0.9090909 | 131 | 97–99 | 0.023359 | 2,620 | 131 | 0.260 |
| 120; 1.0 | 10 | ring | 2 | 0.8941690 | 0.9559229 | 277 | 60–74 | 6.5484e-11 | 6,094 | 277 | 1.340 |
| 120; 1.0 | 10 | circulant4 | 4 | 0.4795755 | 0.7355372 | 41 | 12–13 | 3.638e-10 | 1,804 | 41 | 0.382 |
| 120; 1.0 | 10 | complete | 10 | 0.0000000 | 0.0000000 | 1 | 1–1 | 9.0949e-13 | 110 | 1 | 0.002 |
| 120; 1.0 | 50 | star | 50 | 0.9803922 | 0.9803922 | 825 | 637–641 | 0.01289 | 82,500 | 825 | 1.655 |
| 120; 1.0 | 50 | ring | 2 | 0.9949470 | 0.9979495 | 7952 | 2149–2236 | 5.2969e-09 | 811,104 | 2000 | 9.983 |
| 120; 1.0 | 50 | circulant4 | 4 | 0.8572155 | 0.9384852 | 258 | 69–73 | 7.276e-11 | 52,632 | 258 | 2.456 |
| 120; 1.0 | 50 | complete | 50 | 0.0000000 | 0.0000000 | 1 | 1–1 | 0 | 2,550 | 1 | 0.002 |
| 120; 1.0 | 100 | star | 100 | 0.9900990 | 0.9900990 | 1811 | 1410–1412 | 0.0094128 | 362,200 | 1811 | 3.694 |
| 120; 1.0 | 100 | ring | 2 | 0.9987104 | 0.9994772 | 34450 | 7976–9482 | 3.4925e-10 | 6,958,900 | 2000 | 10.041 |
| 120; 1.0 | 100 | circulant4 | 4 | 0.9242894 | 0.9683168 | 560 | 151–158 | 4.6566e-10 | 226,240 | 560 | 5.338 |
| 120; 1.0 | 100 | complete | 100 | 0.0000000 | 0.0000000 | 1 | 1–1 | 0 | 10,100 | 1 | 0.003 |
| 120; 1.0 | 200 | star | 200 | 0.9950249 | 0.9950249 | 3952 | 3083–3093 | 0.0068848 | 1,580,800 | 2000 | 4.468 |
| 120; 1.0 | 200 | ring | 2 | 0.9996743 | 0.9998680 | 149272 | 40583–41156 | 0 | 60,007,344 | 2000 | 10.457 |
| 120; 1.0 | 200 | circulant4 | 4 | 0.9621072 | 0.9843964 | 1254 | 318–326 | 3.7253e-09 | 1,008,216 | 1254 | 12.189 |
| 120; 1.0 | 200 | complete | 200 | 0.0000000 | 0.0000000 | 1 | 1–1 | 0 | 40,200 | 1 | 0.003 |
| 10000; 0.001 | 10 | star | 10 | 0.9090909 | 0.9090909 | 178 | 145–148 | 0.026817 | 3,560 | 178 | 0.358 |
| 10000; 0.001 | 10 | ring | 2 | 0.8941690 | 0.9559229 | 376 | 110–114 | 8.7311e-10 | 8,272 | 376 | 1.859 |
| 10000; 0.001 | 10 | circulant4 | 4 | 0.4795755 | 0.7355372 | 56 | 18–18 | 2.3283e-10 | 2,464 | 56 | 0.514 |
| 10000; 0.001 | 10 | complete | 10 | 0.0000000 | 0.0000000 | 1 | 1–1 | 0 | 110 | 1 | 0.002 |
| 10000; 0.001 | 50 | star | 50 | 0.9803922 | 0.9803922 | 1048 | 863–867 | 0.013658 | 104,800 | 1048 | 2.146 |
| 10000; 0.001 | 50 | ring | 2 | 0.9949470 | 0.9979495 | 10106 | 2978–3091 | 5.0664e-07 | 1,030,812 | 2000 | 9.856 |
| 10000; 0.001 | 50 | circulant4 | 4 | 0.8572155 | 0.9384852 | 327 | 98–102 | 7.4506e-09 | 66,708 | 327 | 3.087 |
| 10000; 0.001 | 50 | complete | 50 | 0.0000000 | 0.0000000 | 1 | 1–1 | 0 | 2,550 | 1 | 0.002 |
| 10000; 0.001 | 100 | star | 100 | 0.9900990 | 0.9900990 | 2255 | 1851–1853 | 0.0091535 | 451,000 | 2000 | 4.100 |
| 10000; 0.001 | 100 | ring | 2 | 0.9987104 | 0.9994772 | 42906 | 12000–12587 | 3.7253e-08 | 8,667,012 | 2000 | 9.893 |
| 10000; 0.001 | 100 | circulant4 | 4 | 0.9242894 | 0.9683168 | 697 | 206–214 | 4.4703e-08 | 281,588 | 697 | 6.714 |
| 10000; 0.001 | 100 | complete | 100 | 0.0000000 | 0.0000000 | 1 | 1–1 | 7.4506e-09 | 10,100 | 1 | 0.003 |

## Reliable loss and delay model

Each hypothetical data-plus-ack attempt occupies a 1 or 10 ms slot. Directed attempts fail independently with probability 0, .05, or .20. Successful messages are buffered; only failed messages retry. Every node waits for completion of every message in the round before applying the original synchronous update. A 120-round event-level test at 20% loss verifies every state bitwise against the lossless reference. One private field-element exchange per operator-cycle edge is charged before averaging; channel establishment and cryptographic overhead are excluded.

Two schedules expose the hardware assumption: parallel_links permits all links simultaneously; one_peer_per_node uses deterministic greedy edge-color batches, permitting a node to exchange with only one peer at a time. The latter charges the star hub and complete graph for their high degree. It is a specific feasible schedule, not a claim of optimal scheduling.

The network Monte Carlo samples the maximum-of-geometric retry law per batch, with a recorded per-episode tail truncation bound below 1e-12. Expected attempts equal delivered scalars/(1-loss). Deadline fractions are finite-sample estimates, with Wilson 95% intervals in JSON; a fraction of 1 is not a probability-one guarantee. The 1, 10, 60, and 300 second deadlines are illustrative coordination budgets; the additional 900-second deadline matches the separate physical study's 15-minute interval. On a missed deadline this protocol makes no output available: it does not decode early. No physical fallback or response to a missed bid is simulated here. The model excludes correlated/burst loss, disconnects, malicious behavior, finite retry limits, clock synchronization, bandwidth/packet-size effects, encryption, and control traffic. Persistent loss can prevent progress; no deadline guarantee follows from the unbounded retry model.

Selected fine-grid 20%-loss, 1 ms-slot cases below report expected network time and the Monte Carlo 95th percentile. All 0/.05/.20-loss and 1/10 ms-slot results, including 1-second deadline success, are in the JSON. Host kernel times are kept separate and are not added to synthetic network time.

| Operators | Graph | Schedule | Batches/round | Expected seconds | Sample p95 seconds | Fraction <=10 s | Fraction <=60 s | Expected scalar attempts |
|---:|:---|:---|---:|---:|---:|---:|---:|---:|
| 10 | star | parallel_links | 1 | 0.489 | 0.508 | 1.000 | 1.000 | 4,462.5 |
| 10 | star | one_peer_per_node | 10 | 2.600 | 2.650 | 1.000 | 1.000 | 4,462.5 |
| 10 | ring | parallel_links | 1 | 1.052 | 1.078 | 1.000 | 1.000 | 10,352.5 |
| 10 | ring | one_peer_per_node | 3 | 2.300 | 2.347 | 1.000 | 1.000 | 10,352.5 |
| 10 | circulant4 | parallel_links | 1 | 0.182 | 0.193 | 1.000 | 1.000 | 3,092.5 |
| 10 | circulant4 | one_peer_per_node | 5 | 0.630 | 0.652 | 1.000 | 1.000 | 3,092.5 |
| 10 | complete | parallel_links | 1 | 0.006 | 0.008 | 1.000 | 1.000 | 150.0 |
| 10 | complete | one_peer_per_node | 15 | 0.036 | 0.041 | 1.000 | 1.000 | 150.0 |
| 50 | star | parallel_links | 1 | 3.905 | 3.947 | 1.000 | 1.000 | 131,062.5 |
| 50 | star | one_peer_per_node | 50 | 76.422 | 76.687 | 0.000 | 0.000 | 131,062.5 |
| 50 | ring | parallel_links | 1 | 37.746 | 37.887 | 0.000 | 1.000 | 1,288,577.5 |
| 50 | ring | one_peer_per_node | 3 | 81.396 | 81.625 | 0.000 | 0.000 | 1,288,577.5 |
| 50 | circulant4 | parallel_links | 1 | 1.365 | 1.391 | 1.000 | 1.000 | 83,447.5 |
| 50 | circulant4 | one_peer_per_node | 6 | 5.486 | 5.546 | 1.000 | 1.000 | 83,447.5 |
| 50 | complete | parallel_links | 1 | 0.009 | 0.011 | 1.000 | 1.000 | 3,250.0 |
| 50 | complete | one_peer_per_node | 63 | 0.205 | 0.216 | 1.000 | 1.000 | 3,250.0 |
| 100 | star | parallel_links | 1 | 9.364 | 9.429 | 1.000 | 1.000 | 563,875.0 |
| 100 | star | one_peer_per_node | 100 | 328.861 | 329.394 | 0.000 | 0.000 | 563,875.0 |
| 100 | ring | parallel_links | 1 | 178.380 | 178.656 | 0.000 | 0.000 | 10,833,890.0 |
| 100 | ring | one_peer_per_node | 3 | 382.029 | 382.498 | 0.000 | 0.000 | 10,833,890.0 |
| 100 | circulant4 | parallel_links | 1 | 3.203 | 3.239 | 1.000 | 1.000 | 352,110.0 |
| 100 | circulant4 | one_peer_per_node | 6 | 12.635 | 12.719 | 0.000 | 1.000 | 352,110.0 |
| 100 | complete | parallel_links | 1 | 0.010 | 0.012 | 1.000 | 1.000 | 12,750.0 |
| 100 | complete | one_peer_per_node | 127 | 0.461 | 0.477 | 1.000 | 1.000 | 12,750.0 |

## Interpretation and manuscript wording

Topology changes both the spectral convergence rate and per-round communication. A complete graph reaches the average in one synchronous update but requires quadratic total edges and high per-node degree. A ring has bounded degree but slow convergence. The degree-at-most-four circulant adds long-range links and improves convergence without a growing degree; the utility-centered star retains a degree-n bottleneck. These conclusions concern the specified graph families and gain, not arbitrary networks.

Suggested text: 'A reproducible communication study uses 3, 10, 50, 100 and 200 operators and compares star, ring, degree-at-most-four circulant, and complete public consensus graphs. Fixed public iteration budgets combine conservative spectral and binary64 error bounds; all encoded totals are recovered across three fresh public-mask trials per case. Offline decoding thresholds, message counts, and one-host computation times are reported separately. Reliable buffered retransmission preserves the synchronous iteration under independent modeled losses, while changing the completion time. A node-limited matching schedule exposes the star hub degree and distinguishes connectivity from practical communication scalability. Slot duration and deadlines are scenario assumptions, not measured network or microgrid timings.'

The public communication graph and private mask graph remain distinct. Improving the former does not improve the latter's collusion tolerance. This study does not extend privacy to constrained bids, metering, active attacks, or dropped participants.

Artifacts: `communication_validation.json`, `convergence_n100.svg`, and this report. JSON contains complete inputs, masks, seeds, bounds, sample traces, network assumptions, hardware/software metadata and all results.
