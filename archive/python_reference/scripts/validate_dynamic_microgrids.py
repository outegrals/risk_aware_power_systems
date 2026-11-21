#!/usr/bin/env python3
"""Synthetic sequential microgrid evidence, not an electrical transient model.

Run: python3 scripts/validate_dynamic_microgrids.py
Only NumPy plus the standard library are required. No downloaded measurements.
Physical controls never read bids. Conditional moments are evaluated before
drawing the next realization, exactly over a declared finite distribution.
"""
from __future__ import annotations
import csv
import hashlib
import json
import platform
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/microgrid"
DT, STEPS, REPS = .25, 192, 30
SCALES, SIZES = (0., .5, 1., 2.), (10, 50, 100)
SEED = 2026111800
Z, W = np.array([-1., 0., 1.]), np.array([.2, .6, .2])
JOINT = (W[:, None]*W[None, :]).reshape(9)
ETA, DEGRADATION, A1 = .95, 5., 20.
NAMES = ("RN", "fixed_ERS", "retuned")


def assets(n, seed):
    rng = np.random.default_rng(seed)
    size = np.exp(rng.uniform(np.log(.35), np.log(2.), n))
    return {"size": size, "capacity": 4*size, "power": .65*size,
            "export_limit": 1.8*size, "phase": rng.uniform(-.6, .6, n)}


def conditional_physics(a, energy, step, scale):
    """Shape (common state, independent local state, operator), before sampling."""
    hour = (step*DT) % 24
    size, phase = a["size"], a["phase"]
    wind = 1.8*size*(.62+.10*np.sin(2*np.pi*hour/24+phase))
    solar = 2.4*size*max(0., np.sin(np.pi*(hour-6)/12))
    local_load = .65*size*(.85+.08*np.cos(2*np.pi*hour/24+phase)
                          +.18*np.exp(-((hour-19)/3)**2))
    common, local = Z[:, None, None], Z[None, :, None]
    gen = (wind+solar)*(1+scale*(.18*common+.12*local))
    load = local_load*(1+scale*(.04*common+.04*local))
    surplus = gen-load
    reference = size*(.60+.10*np.sin(2*np.pi*hour/24))
    max_charge = np.minimum(a["power"], (.9*a["capacity"]-energy)/(ETA*DT))
    max_discharge = np.minimum(a["power"], ETA*(energy-.1*a["capacity"])/DT)
    charge = np.minimum(np.maximum(surplus-reference, 0), np.maximum(0, max_charge))
    discharge = np.minimum(np.maximum(reference-surplus, 0), np.maximum(0, max_discharge))
    net = surplus+discharge-charge
    exports = np.clip(net, 0, a["export_limit"])
    imports, spill = np.maximum(-net, 0), np.maximum(net-a["export_limit"], 0)
    next_energy = energy+DT*(ETA*charge-discharge/ETA)
    external_base = 2.7*size.sum()*(1+.10*np.cos(2*np.pi*(hour-19)/24))
    external_cond_mean = external_base*(1+.06*scale*Z)
    return {"generation": gen, "load": load, "export": exports,
            "import": imports, "spill": spill, "charge": charge,
            "discharge": discharge, "next_energy": next_energy,
            "external_base": external_base, "external_cond_mean": external_cond_mean}


def mean(values):
    return np.einsum("a,b,abi->i", W, W, values)


def inputs_and_constant(p, a2):
    x, imp = p["export"], p["import"]
    mu = mean(x)
    # Lower 0.1 quantile: F(u-) <= 0.1, then cap at mean and feeder limit.
    flat = x.reshape(9, -1)
    order = np.argsort(flat, axis=0, kind="stable")
    sorted_values = np.take_along_axis(flat, order, axis=0)
    weights = JOINT[order]
    index = np.argmax(np.cumsum(weights, axis=0) >= .1, axis=0)
    caps = np.minimum(mu, sorted_values[index, np.arange(len(mu))])
    signed = x-imp
    cond_mean = np.einsum("b,abi->ai", W, signed)
    cond_var = np.einsum("b,abi->ai", W, (signed-cond_mean[:, None, :])**2)
    residual_mean = p["external_cond_mean"]-cond_mean.sum(axis=1)
    # E[(a1+2*a2*(external-sum signed))*signed_i], including own covariance.
    expected_rt_net = A1*mean(signed) + 2*a2*np.einsum(
        "a,ai->i", W, residual_mean[:, None]*cond_mean-cond_var)
    constant = expected_rt_net-DEGRADATION*mean(p["charge"]+p["discharge"])
    forecast_load = p["external_base"]+mean(imp).sum()
    shortfall_probability = np.einsum("a,b,abi->i", W, W, x < caps[None, None, :]-1e-12)
    assert shortfall_probability.max() <= .1+1e-12
    return mu, caps, constant, forecast_load


def level_solution(total, caps, inverse_k):
    """Solve inverse_k*t + sum min(cap_i,t) = total by sorted breakpoints."""
    ordered = np.sort(caps)
    prefix = np.r_[0., np.cumsum(ordered)]
    n = len(caps)
    if inverse_k == 0 and total >= prefix[-1]-1e-12:
        return ordered[-1]
    for j in range(n+1):
        denominator = inverse_k+n-j
        if denominator <= 0:
            continue
        level = max(0., (total-prefix[j])/denominator)
        if j == n or level <= ordered[j]+1e-12:
            return level
    raise AssertionError("no scalar root")


def strategies(mu, caps):
    n, total = len(mu), mu.sum()
    qrn = np.minimum(caps, level_solution(total, caps, 1))
    qers = np.minimum(caps, level_solution(total, caps, n))
    target = min(total/2, caps.sum())
    level = level_solution(target, caps, 0)
    qret = np.minimum(caps, level)
    kret = level/(total-target)
    q = np.stack([qrn, qers, qret])
    k = np.array([1., 1/n, kret])
    derivative = k[:, None]*(total-q.sum(axis=1))[:, None]-q
    br = np.clip(k[:, None]/(1+k[:, None])*(total-q.sum(axis=1)[:, None]+q), 0, caps)
    residual = np.max(np.abs(br-q))/max(1., total)
    assert residual < 1e-10
    assert abs(qret.sum()-target) < 1e-10*max(1., total)
    assert np.all(q >= -1e-12) and np.all(q <= caps+1e-12)
    return q, k, residual


def one_run(n, scale, rep, trace=False):
    a = assets(n, SEED+rep)
    rng = np.random.default_rng(SEED+10000+rep)
    energy = .5*a["capacity"]
    a2 = 5/n
    expected_profit, realized_profit, shortfall, expected_shortfall = (np.zeros(3) for _ in range(4))
    schedule, cap_hits = np.zeros(3), np.zeros(3)
    checks = {"max_normalized_br_residual": 0., "max_balance_residual_MW": 0.,
              "max_soc_violation_MWh": 0., "min_conventional_generation_MW": float("inf"),
              "max_simultaneous_charge_discharge_MW2": 0., "max_power_limit_violation_MW": 0.}
    soc_low, soc_high, imports, curtailment, throughput = 1., 0., 0., 0., 0.
    rows = []
    reference_timings = []
    for step in range(STEPS):
        p = conditional_physics(a, energy, step, scale)
        mu, caps, const, load_mean = inputs_and_constant(p, a2)
        assert mu.min() > 0, "the paper's positive-forecast condition must hold"
        start = time.perf_counter_ns()
        q, k, residual = strategies(mu, caps)
        reference_timings.append((time.perf_counter_ns()-start)/1e6)
        expected = DT*(const.sum()+2*a2*q.sum(axis=1)*(mu.sum()-q.sum(axis=1)))
        expected_profit += expected
        schedule += DT*q.sum(axis=1)
        cap_hits += np.sum(np.abs(q-caps)<1e-10, axis=1)
        expected_shortfall += DT*np.einsum("a,b,sabi->s", W, W,
            np.maximum(q[:, None, None, :]-p["export"][None, :, :, :], 0))
        # Verify the entire conditional support, not only the drawn realization.
        balance = p["generation"]+p["discharge"]+p["import"]-p["load"]-p["charge"]-p["export"]-p["spill"]
        checks["max_balance_residual_MW"] = max(checks["max_balance_residual_MW"], float(abs(balance).max()))
        ev = max(0., float((.1*a["capacity"]-p["next_energy"]).max()),
                 float((p["next_energy"]-.9*a["capacity"]).max()))
        checks["max_soc_violation_MWh"] = max(checks["max_soc_violation_MWh"], ev)
        checks["max_normalized_br_residual"] = max(checks["max_normalized_br_residual"], residual)
        checks["max_simultaneous_charge_discharge_MW2"] = max(checks["max_simultaneous_charge_discharge_MW2"], float((p["charge"]*p["discharge"]).max()))
        checks["max_power_limit_violation_MW"] = max(checks["max_power_limit_violation_MW"], float(np.maximum(0, p["charge"]+p["discharge"]-a["power"]).max()))
        max_signed = (p["export"]-p["import"]).max(axis=1).sum(axis=1)
        min_pg = min(float((p["external_cond_mean"]-.02*scale*p["external_base"]-max_signed).min()), float(load_mean-q.sum(axis=1).max()))
        checks["min_conventional_generation_MW"] = min(checks["min_conventional_generation_MW"], min_pg)
        assert ev < 1e-10 and min_pg > 0
        common = rng.choice(3, p=W)
        local = rng.choice(3, size=n, p=W)
        external_state = rng.choice(3, p=W)
        draw = lambda key: p[key][common, local, np.arange(n)]
        x, imp = draw("export"), draw("import")
        c, d, spill = draw("charge"), draw("discharge"), draw("spill")
        ext = p["external_cond_mean"][common]+.02*scale*p["external_base"]*Z[external_state]
        lambda_rt = A1+2*a2*(ext+imp.sum()-x.sum())
        lambda_da = A1+2*a2*(load_mean-q.sum(axis=1))
        realized_profit += DT*(lambda_da*q.sum(axis=1)+lambda_rt*(x.sum()-imp.sum()-q.sum(axis=1))-DEGRADATION*(c+d).sum())
        shortfall += DT*np.maximum(q-x, 0).sum(axis=1)
        imports += DT*imp.sum()
        curtailment += DT*spill.sum()
        throughput += DT*(c+d).sum()
        energy = draw("next_energy")
        soc_low, soc_high = min(soc_low, float((energy/a["capacity"]).min())), max(soc_high, float((energy/a["capacity"]).max()))
        if trace:
            rows.append({"hour_end": (step+1)*DT, "mean_soc": float(np.mean(energy/a["capacity"])),
                         "min_soc": float(min(energy/a["capacity"])), "max_soc": float(max(energy/a["capacity"])),
                         "mean_export_MW": float(x.mean()), "mean_import_MW": float(imp.mean()),
                         "mean_forecast_MW": float(mu.mean()), "mean_cap_MW": float(caps.mean()),
                         "RN_mean_bid_MW": float(q[0].mean()), "ERS_mean_bid_MW": float(q[1].mean()),
                         "retuned_mean_bid_MW": float(q[2].mean()), "mean_next_energy_MWh": float(energy.mean())})
    return {"n": n, "uncertainty_scale": scale, "replicate": rep, "asset_seed": SEED+rep,
            "shock_seed": SEED+10000+rep, "expected_cash": expected_profit.tolist(),
            "realized_cash": realized_profit.tolist(), "scheduled_MWh": schedule.tolist(),
            "shortfall_MWh": shortfall.tolist(), "expected_shortfall_MWh": expected_shortfall.tolist(),
            "cap_hit_fraction": (cap_hits/(n*STEPS)).tolist(), "imports_MWh": imports,
            "curtailment_MWh": curtailment, "throughput_MWh": throughput,
            "soc_range": [soc_low, soc_high], "initial_energy_MWh": float((.5*a["capacity"]).sum()),
            "terminal_energy_MWh": float(energy.sum()), "validation": checks,
            "reference_solver_ms_median": float(np.median(reference_timings)),
            "reference_solver_ms_p95": float(np.quantile(reference_timings, .95))}, rows


def aggregate(runs):
    result = []
    for n in SIZES:
        for scale in SCALES:
            selected = [r for r in runs if r["n"] == n and r["uncertainty_scale"] == scale]
            cash = np.array([r["realized_cash"] for r in selected])
            gain = (cash[:, 2]-cash[:, 0])/n
            # Replicates, not intervals, are the independent sampling units.
            mean_gain = float(gain.mean())
            margin = float(2.04523*gain.std(ddof=1)/np.sqrt(REPS))
            result.append({"n": n, "scale": scale,
                "cash_per_operator_48h": cash.mean(axis=0).tolist() if n == 1 else (cash.mean(axis=0)/n).tolist(),
                "expected_cash_per_operator_48h": (np.array([r["expected_cash"] for r in selected]).mean(axis=0)/n).tolist(),
                "retuned_minus_RN_cash_per_operator": mean_gain,
                "paired_95pct_t_interval": [mean_gain-margin, mean_gain+margin],
                "shortfall_kWh_per_operator_48h": (1000*np.array([r["shortfall_MWh"] for r in selected]).mean(axis=0)/n).tolist(),
                "expected_shortfall_kWh_per_operator_48h": (1000*np.array([r["expected_shortfall_MWh"] for r in selected]).mean(axis=0)/n).tolist(),
                "cap_hit_fraction": np.array([r["cap_hit_fraction"] for r in selected]).mean(axis=0).tolist(),
                "mean_imports_kWh_per_operator": 1000*float(np.mean([r["imports_MWh"] for r in selected]))/n,
                "mean_curtailment_kWh_per_operator": 1000*float(np.mean([r["curtailment_MWh"] for r in selected]))/n,
                "mean_terminal_minus_initial_kWh_per_operator": 1000*float(np.mean([r["terminal_energy_MWh"]-r["initial_energy_MWh"] for r in selected]))/n})
    return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    runs, trace = [], []
    for n in SIZES:
        for scale in SCALES:
            for rep in range(REPS):
                result, rows = one_run(n, scale, rep, trace=(n == 50 and scale == 1 and rep == 0))
                runs.append(result)
                trace.extend(rows)
            print(f"completed n={n} scale={scale}", flush=True)
    summary = aggregate(runs)
    report = {"design": {"synthetic": True, "dt_hours": DT, "steps": STEPS, "replicates": REPS,
              "sizes": SIZES, "scales": SCALES, "seed": SEED, "strategy_order": NAMES,
              "conditional_state_values": Z.tolist(), "conditional_probabilities": W.tolist(),
              "a1": A1, "a2": "5/n", "efficiency_each_direction": ETA,
              "throughput_cost_per_MWh": DEGRADATION, "cap": "min(mean export, lower 0.1 quantile export)",
              "forecast_horizon": "one 15-minute interval; not a fixed prior-day 24-hour schedule",
              "uncertainty": "common and local discrete shocks independent over time; SOC carries history",
              "limitations": "Ideal current inputs; cap solver uses global information; no constrained privacy claim; no voltage/frequency/network-flow model; no terminal SOC target or battery optimization."},
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "environment": {"python": platform.python_version(), "numpy": np.__version__,
                              "platform": platform.platform(), "machine": platform.machine()},
              "elapsed_seconds": time.perf_counter()-start, "summary": summary, "runs": runs}
    (OUT/"results.json").write_text(json.dumps(report, indent=2)+"\n")
    with (OUT/"representative_trajectory.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(trace[0]))
        writer.writeheader()
        writer.writerows(trace)
    lines = ["# Synthetic microgrid study", "", "30 paired 48-hour trajectories per setting; 15-minute steps. Cash in synthetic currency units; no real data calibration. Each confidence interval is an unadjusted, approximate paired Student-t interval across independent replicates, not a simultaneous familywise claim.", "", "| n | s | RN cash/operator | fixed ERS | retuned | gain [95% interval] | shortfall kWh/operator RN / ERS / retuned |", "|---:|---:|---:|---:|---:|---|---|"]
    for row in summary:
        c = row["cash_per_operator_48h"]
        lo, hi = row["paired_95pct_t_interval"]
        sf = row["shortfall_kWh_per_operator_48h"]
        lines.append(f"| {row['n']} | {row['scale']:g} | {c[0]:.2f} | {c[1]:.2f} | {c[2]:.2f} | {row['retuned_minus_RN_cash_per_operator']:.2f} [{lo:.2f}, {hi:.2f}] | {sf[0]:.2f} / {sf[1]:.2f} / {sf[2]:.2f} |")
    lines += ["", "All strategies share the same bid-independent physical controller and therefore the same SOC, imports, curtailment, and battery throughput. Shortfall measures contracted export missing at settlement, not unserved local load. Import purchases and battery throughput costs are included in cash. No terminal energy value or recharge requirement is imposed; identical terminal states make the paired comparison valid but absolute cash is not optimized lifetime profit.", "", f"Wall time for this validation: {report['elapsed_seconds']:.2f} s; this includes validation and physical simulation and is not a distributed deployment benchmark."]
    (OUT/"validation_report.md").write_text("\n".join(lines)+"\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
