# Synthetic microgrid energy-management validation

This directory contains a sequential, grid-connected microgrid application of the paper's bounded bidding game. It models battery energy states and interval energy accounting, not voltage/frequency transients, converter controls, network power flow, islanding, or an optimized battery policy. All profiles and parameters are synthetic; no measured data or market calibration is used.

## Reproduce and inspect

From the repository root, with Python 3 and NumPy available:

```sh
python3 scripts/validate_dynamic_microgrids.py
```

The script writes:

- `output/microgrid/results.json`: parameters, source SHA-256, Python/NumPy/platform versions, all 360 runs, aggregate statistics, and physical/equilibrium checks.
- `output/microgrid/validation_report.md`: the comparison table and interpretation limits.
- `output/microgrid/representative_trajectory.csv`: interval summaries for `n=50`, uncertainty scale `1`, replicate `0`. This is one representative path, not the full ensemble.

The numerical outputs are reproducible from the recorded seeds. Timing fields depend on the machine and execution conditions; they measure the centralized reference calculation and simulation, not a deployed distributed controller. Rerunning overwrites these three generated files; this README is maintained separately.

## Experimental design and statistical unit

| Parameter | Implemented value |
|---|---|
| Microgrid count `n` | 10, 50, 100 |
| Raw disturbance amplitude `s` | 0, 0.5, 1, 2 |
| Independent replicates per setting | 30 |
| Horizon | 192 intervals, 48 hours |
| Interval duration `dt` | 0.25 hours (15 minutes) |
| Base seed | 2026111800 |
| Asset seed, replicate `r` | base seed + `r` |
| Shock seed, replicate `r` | base seed + 10000 + `r` |
| Strategy order in JSON arrays | RN, fixed_ERS, retuned |

The 15-minute interval is a **synthetic supervisory study choice**, not an observed market cadence or a simulation of a fixed schedule submitted the preceding day. At each interval, the current battery state is observed and the next interval's finite conditional distribution is evaluated before sampling its realization. The forecasts are ideal conditional model means; forecast-model estimation error is not tested.

Within each replicate, all strategies share assets, shocks, and the physical controller. Thus strategy differences are paired. Replicates, not individual intervals or microgrids, are the independent statistical units. Asset draws as well as trajectory shocks vary between replicates; even `s=0` comparisons vary because assets differ. Seeds are reused across settings, so different uncertainty/population settings are not asserted to be independent or paired physical systems.

The reported gain interval is `mean(gain) +/- 2.04523 * sd(gain)/sqrt(30)`, where gain is the 48-hour retuned-minus-RN cash difference per operator. It is an approximate, unadjusted 95% paired Student-t interval for each setting, not a simultaneous familywise guarantee. No normality of individual interval cash is assumed or established. The sum of exact one-interval conditional expected cash along each sampled SOC path is also recorded, separately from realized trajectory cash; it is not an exact unconditional 48-hour expectation.

## Assets, profiles, and uncertainty

For microgrid `i`, the size factor `b_i` is log-uniform on `[0.35, 2]`, and phase `phi_i` is uniform on `[-0.6, 0.6]` radians. The script draws all sizes and then all phases with the asset generator. Capacities scale with `b_i`:

| Asset quantity | Value |
|---|---|
| Battery energy capacity `C_i` | `4 b_i` MWh |
| Charge and discharge power limit | `0.65 b_i` MW each |
| Export connection limit | `1.8 b_i` MW |
| Initial energy | `0.5 C_i` |
| Energy bounds | `0.1 C_i <= E_i <= 0.9 C_i` |
| Charging/discharging efficiency | `eta_c = eta_d = 0.95` |
| Synthetic throughput cost | 5 cash units/MWh of `charge + discharge` |

With clock hour `h=(step*dt) mod 24`, nominal profiles in MW are:

```text
wind_i(h)  = 1.8 b_i [0.62 + 0.10 sin(2 pi h/24 + phi_i)]
solar_i(h) = 2.4 b_i max(0, sin(pi (h-6)/12))
load_i(h)  = 0.65 b_i [0.85 + 0.08 cos(2 pi h/24 + phi_i)
                     + 0.18 exp(-((h-19)/3)^2)]
```

Each independent shock takes values `{-1,0,1}` with probabilities `{0.2,0.6,0.2}`. There is one common shock `z_c` shared by all microgrids, one independent local shock `z_i` for each microgrid, and one independent external-load shock `z_e`. At scale `s`:

```text
generation_i = (wind_i + solar_i) [1 + s(0.18 z_c + 0.12 z_i)]
local_demand_i = load_i [1 + s(0.04 z_c + 0.04 z_i)]
external_base = 2.7 sum_i b_i [1 + 0.10 cos(2 pi (h-19)/24)]
external_demand = external_base [1 + s(0.06 z_c + 0.02 z_e)]
```

Common shocks create cross-microgrid dependence; shared common/local shocks also correlate each microgrid's generation and demand. Conditional on the common shock, local shocks are independent across microgrids. Shocks are independent over time; battery energy carries history. This is not a temporally correlated weather model.

Each operator's conditional export distribution is computed over its nine common/local outcomes with product probabilities. Aggregate conditional cash moments account for the common dependence and each operator's own covariance; the implementation does not incorrectly treat operators as unconditionally independent. The independent external shock has zero conditional mean. Generation and local-demand multipliers stay positive at every tested scale. Storage limits and export clipping can change conditional means and caps as `s` changes, so this experiment does **not** hold means and feasible sets fixed as in the separate covariance-only sensitivity test.

## Bid-independent physical controller and accounting

The physical policy tracks the synthetic reference

```text
reference_i(h) = b_i [0.60 + 0.10 sin(2 pi h/24)].
```

It charges from surplus above this reference, or discharges to cover a deficit below it. Its permissible interval powers are

```text
charge_max_i = min(0.65 b_i, (0.9 C_i-E_i)/(eta_c dt))
discharge_max_i = min(0.65 b_i, eta_d (E_i-0.1 C_i)/dt).
```

The realized charge/discharge command is the corresponding nonnegative reference error capped by the nonnegative permissible power. The commands cannot be positive simultaneously. With surplus equal to generation minus local demand:

```text
net_i = surplus_i + discharge_i - charge_i
export_i = clip(net_i, 0, 1.8 b_i)
import_i = max(-net_i, 0)
curtailment_i = max(net_i-1.8 b_i, 0)
E_i_next = E_i + dt (eta_c charge_i - discharge_i/eta_d).
```

Grid imports supply deficits without a modeled import limit. The export limit is not a general feeder thermal/voltage model. The checked local balance is

```text
generation + discharge + import
    = local demand + charge + export + curtailment.
```

External demand excludes the local demands already subtracted above. Consequently conventional generation is `external_demand + sum(import) - sum(export)`, avoiding double counting. The script checks positive conventional generation over the complete next-interval support and the pre-interval scheduling balance. It records maximum local-balance, SOC, simultaneous-charge/discharge, and power-limit residuals over the conditional support, not just sampled realizations.

The physical controller never reads bids or strategy coefficients. All strategies therefore have identical SOC, imports, curtailment, and battery throughput on a paired path. This preserves bid-exogenous physical output for the interval game; it does not demonstrate optimal storage operation. The scenarios are predominantly exporting portfolios, so they provide limited evidence about substantial importing operation.

No terminal SOC target, replenishment obligation, or terminal energy value is imposed. Identical terminal states prevent a differential battery-depletion advantage in the paired strategy comparison, but absolute cash is not optimized lifetime profit. Initial and terminal energy, throughput, imports, and curtailment are recorded explicitly.

## Conditional forecasts, caps, controllers, and settlement

Let `mu_i = E[export_i | current information]`, `X = sum(mu_i)`, and `Q = sum(q_i)`. All conditional means are computed before drawing the new shocks. Forecasts are checked to remain strictly positive. The frozen bidding cap is

```text
u_i = min(mu_i, lower 0.1 quantile of conditional export_i).
```

The quantile is the first sorted value whose cumulative probability reaches 0.1. Export outcomes already obey the connection limit. The script verifies `Pr(export_i < u_i) <= 0.1`, with a `1e-12` comparison tolerance. Since every bid is at most its cap, this is an **individual, one-interval conditional export-shortfall bound** under the known finite law. It is not a simultaneous all-operator/all-interval reliability guarantee, a guarantee of realized delivery, or a bound on shortfall magnitude.

With `k_i = 1 + 4 a2 beta_i`, the strategies are:

- **RN:** fixed common `k=1`.
- **Fixed ERS:** the original uncapped common coefficient `k=1/n`, solving the actual capped game at that coefficient.
- **Retuned:** target `Q_target=min(X/2,sum(u))`; find `t` with `sum(min(u_i,t))=Q_target`, then use common `k=t/(X-Q_target)`. When all caps bind, the implementation chooses `t=max(u)`.

Fixed-coefficient equilibria solve `r + sum(min(u_i,k r)) = X`, `q_i=min(u_i,k r)`. The implementation uses sorted breakpoints and checks the result against each direct unilateral quadratic best response, bid bounds, and the retuned aggregate target. The cap solver has global input access. Retuning is a coordinated controller design; it is not a proved equilibrium in strategically chosen `beta`.

Synthetic market coefficients are `a1=20` and `a2=5/n`. Thus population comparisons scale the external load and physical assets and normalize the price slope; they are not tests of an unchanged physical market. Prices and cash are uncalibrated synthetic currency quantities:

```text
forecast_market_load = external_base + sum(E[import_i])
price_pre = a1 + 2 a2 (forecast_market_load - Q)
price_real = a1 + 2 a2 (external_demand + sum(import_i) - sum(export_i))
cash = dt [price_pre Q
           + price_real (sum(export_i) - sum(import_i) - Q)
           - 5 sum(charge_i + discharge_i)].
```

The import purchase is included at the realized price; the throughput cost is charged on both charging and discharging at the AC-side powers. No additional capital, replacement, startup, or emissions cost is modeled. Conditional expected cash includes the applicable export/import covariance terms. Multiplication by `dt` converts interval power revenue rates to cash.

Reported shortfall energy is `dt * sum(max(q_i-export_i,0))`: missing contracted exports, not unserved local load. Imports are reported separately. Fixed and retuned ERS can differ in shortfall because they submit different quantities; the study does not establish dominance on every risk measure.

## Why communication evidence is separate

`scripts/validate_communication.py` studies the masked protocol for quantized aggregate forecasts and the **uncapped** ERS allocation. Its communication graphs, arithmetic certificate, private-mask assumptions, and modeled deadlines concern that specific protocol. This microgrid script instead computes conditional quantile caps and a constrained allocation with a centralized reference solver; it neither executes the communication protocol nor charges its deadlines to these physical trajectories.

The two modules therefore cannot be combined into a claim of a private, distributed, deadline-certified constrained microgrid controller. Caps, saturation, and constrained allocation may require extra disclosures not covered by the uncapped privacy theorem. Quantization and feasible-cap effects also need joint treatment. A missed aggregation deadline and its executed physical fallback have not been integrated into this battery study. The separate modules support separate, explicitly scoped claims.

## Modeling sources and evidence limits

These primary sources motivate the control-layer distinction and battery modeling ingredients; the numerical parameters above are the script's synthetic assumptions, not parameters calibrated from these sources.

- Josep M. Guerrero, Juan C. Vasquez, José Matas, Luis García de Vicuña, and Miguel Castilla, **“Hierarchical Control of Droop-Controlled AC and DC Microgrids—A General Approach Toward Standardization,”** *IEEE Transactions on Industrial Electronics*, 58(1), 158–172, 2011. DOI: [10.1109/TIE.2010.2066534](https://doi.org/10.1109/TIE.2010.2066534). The primary/secondary/tertiary distinction supports treating market coordination as a supervisory application, without claiming electrical stability. [IEEE record](https://ieeexplore.ieee.org/document/5546958).
- Nicholas DiOrio, Aron Dobos, Steven Janzou, Austin Nelson, and Blake Lundstrom, **“Technoeconomic Modeling of Battery Energy Storage in SAM,”** NREL/TP-6A20-64641, September 2015. DOI: [10.2172/1225314](https://doi.org/10.2172/1225314). The report documents separate charge/discharge efficiencies, energy-state and rate limits, and dispatch accounting. The present script uses a simpler constant-efficiency model, not the full SAM voltage/thermal/aging model. [OSTI record](https://www.osti.gov/biblio/1225314), [NREL report](https://docs.nrel.gov/docs/fy15osti/64641.pdf).

This is finite synthetic numerical evidence for conditional interval bidding with an evolving physical energy state. It does not establish performance on measured microgrids, multi-period strategic optimality, forecast robustness under misspecification, controller commitment, electrical network feasibility, or end-to-end private distributed deployment.

Independent verification: `python3 scripts/audit_dynamic_microgrids.py` writes `independent_audit.json`. Generate the vector trajectory figure with `python3 scripts/render_microgrid_evidence.py`; the expanded manuscript includes the resulting `trajectory_plot.tex`.
