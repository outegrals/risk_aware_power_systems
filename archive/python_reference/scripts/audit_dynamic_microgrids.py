#!/usr/bin/env python3
"""Independent small-system audit of the sequential microgrid study.

This does not repeat the trajectory experiment. It enumerates the complete
81-state joint law for two operators, evaluates settlements directly, fits
unilateral objectives from those evaluations, and compares the study's sorted
root solver with an independent bisection implementation. Only this script's
own JSON artifact is written.
"""
from __future__ import annotations

import hashlib
import importlib.util
from itertools import product
import json
from pathlib import Path
import sys

import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/validate_dynamic_microgrids.py"
OUTPUT = ROOT / "output/microgrid/independent_audit.json"
PROBABILITIES = np.array([.2, .6, .2])
STATES = np.array([-1., 0., 1.])


def load_source():
    spec = importlib.util.spec_from_file_location("microgrid_study_under_audit", SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def enumerate_joint_law(study, physical, scale):
    """Enumerate common, both independent local, and external-load shocks."""
    probabilities, exports, imports, costs, external_load = [], [], [], [], []
    for common, local0, local1, external in product(range(3), repeat=4):
        probability = np.prod(PROBABILITIES[[common, local0, local1, external]])
        locals_ = [local0, local1]
        x = np.array([physical["export"][common, locals_[i], i] for i in range(2)])
        imp = np.array([physical["import"][common, locals_[i], i] for i in range(2)])
        throughput = np.array([
            physical["charge"][common, locals_[i], i]
            + physical["discharge"][common, locals_[i], i] for i in range(2)
        ])
        ext = (physical["external_cond_mean"][common]
               + .02 * scale * physical["external_base"] * STATES[external])
        probabilities.append(probability)
        exports.append(x)
        imports.append(imp)
        costs.append(study.DEGRADATION * throughput)
        external_load.append(ext)
    weights = np.array(probabilities)
    assert abs(weights.sum()-1) < 1e-14
    return weights, np.array(exports), np.array(imports), np.array(costs), np.array(external_load)


def settlement_evaluator(study, law, a2):
    weights, exports, imports, costs, external = law
    load = external + imports.sum(axis=1)
    forecast_load = float(weights @ load)
    lambda_rt = study.A1 + 2*a2*(load-exports.sum(axis=1))

    def evaluate(q, beta):
        q = np.array(q)
        lambda_da = study.A1 + 2*a2*(forecast_load-q.sum())
        # Both settlement cash and the price-deviation objective receive DT.
        scenario_cash = study.DT * (
            lambda_da*q + lambda_rt[:, None]*(exports-imports-q)-costs
        )
        expected_cash = weights @ scenario_cash
        squared_spread = float(weights @ ((lambda_rt-lambda_da)**2))
        objective = expected_cash-study.DT*np.array(beta)*squared_spread
        return expected_cash, objective, squared_spread

    return evaluate, forecast_load


def bisection_equilibrium(total, caps, k):
    left, right = 0., float(total)
    for _ in range(100):
        mid = (left+right)/2
        residual = mid+np.minimum(caps, k*mid).sum()-total
        if residual > 0:
            right = mid
        else:
            left = mid
    return np.minimum(caps, k*(left+right)/2)


def audit_state(study, asset_seed, fractions, step, scale):
    assets = study.assets(2, asset_seed)
    energy = assets["capacity"]*np.array(fractions)
    p = study.conditional_physics(assets, energy, step, scale)
    a2 = 5/2
    mu, caps, constant, forecast_load = study.inputs_and_constant(p, a2)
    law = enumerate_joint_law(study, p, scale)
    weights, exports, imports, _, _ = law
    evaluate, independent_load = settlement_evaluator(study, law, a2)
    independent_mu = weights @ exports
    # Recompute each marginal quantile from the complete joint law.
    independent_caps = []
    for i in range(2):
        ordering = np.argsort(exports[:, i], kind="stable")
        cdf = np.cumsum(weights[ordering])
        quantile = exports[ordering[np.flatnonzero(cdf >= .1)[0]], i]
        independent_caps.append(min(independent_mu[i], quantile))
    independent_caps = np.array(independent_caps)
    q, k, _ = study.strategies(mu, caps)
    cash_error, br_error, root_error = 0., 0., 0.
    for profile in [np.zeros(2), caps, *q]:
        cash, _, _ = evaluate(profile, np.zeros(2))
        predicted = study.DT*(constant+2*a2*profile*(mu.sum()-profile.sum()))
        cash_error = max(cash_error, float(abs(cash-predicted).max()))
    for profile, ki in zip(q, k):
        beta = np.full(2, (ki-1)/(4*a2))
        root_error = max(root_error, float(abs(
            profile-bisection_equilibrium(mu.sum(), caps, ki)).max()))
        for i in range(2):
            cap = caps[i]
            if cap < 1e-10:
                direct_br = 0.
            else:
                values = []
                for own_bid in [0., cap/2, cap]:
                    candidate = profile.copy()
                    candidate[i] = own_bid
                    values.append(evaluate(candidate, beta)[1][i])
                j0, jhalf, jcap = values
                quadratic = 2*(jcap-2*jhalf+j0)/cap**2
                linear = (jcap-j0-quadratic*cap**2)/cap
                assert quadratic < 0, (step, scale, fractions, quadratic)
                direct_br = np.clip(-linear/(2*quadratic), 0., cap)
            br_error = max(br_error, abs(float(direct_br-profile[i])))
    shortfall_probabilities = weights @ (exports < caps-1e-12)
    balance = (p["generation"]+p["discharge"]+p["import"]-p["load"]
               -p["charge"]-p["export"]-p["spill"])
    soc_recurrence = (p["next_energy"]-energy
                      -study.DT*(study.ETA*p["charge"]-p["discharge"]/study.ETA))
    result = {
        "asset_seed": asset_seed, "initial_soc_fractions": list(fractions),
        "step": step, "scale": scale, "joint_scenarios": len(weights),
        "conditional_export_means": mu.tolist(), "caps": caps.tolist(),
        "max_mean_error": float(abs(mu-independent_mu).max()),
        "forecast_load_error": abs(float(forecast_load-independent_load)),
        "max_quantile_cap_error": float(abs(caps-independent_caps).max()),
        "max_direct_settlement_cash_error": cash_error,
        "max_direct_objective_best_response_error_MW": br_error,
        "max_bisection_equilibrium_error_MW": root_error,
        "max_cap_shortfall_probability": float(shortfall_probabilities.max()),
        "max_balance_error_MW": float(abs(balance).max()),
        "max_soc_recurrence_error_MWh": float(abs(soc_recurrence).max()),
        "maximum_import_MW": float(imports.max()),
    }
    for name in ["max_mean_error", "forecast_load_error", "max_quantile_cap_error",
                 "max_direct_settlement_cash_error", "max_direct_objective_best_response_error_MW",
                 "max_bisection_equilibrium_error_MW", "max_balance_error_MW",
                 "max_soc_recurrence_error_MWh"]:
        assert result[name] < 1e-9, (name, result)
    assert result["max_cap_shortfall_probability"] <= .1+1e-12, result
    return result


def audit_solver_edges(study):
    cases = [
        ([120., 100., 80.], [30., 20., 10.], "capacity_shortage"),
        ([120., 100., 80.], [50., 50., 50.], "capacity_equals_half_forecast"),
        ([120., 100., 80.], [0., 0., 0.], "all_zero_caps"),
        ([120., 100., 80.], [0., 200., 0.], "one_positive_cap"),
        ([1., 18., 18., 18., 45.], [1., 18., 18., 18., 45.], "heterogeneous_five_operator"),
    ]
    results = []
    for means, limits, name in cases:
        mu, caps = np.array(means), np.array(limits)
        q, k, _ = study.strategies(mu, caps)
        error = max(float(abs(profile-bisection_equilibrium(mu.sum(), caps, ki)).max())
                    for profile, ki in zip(q, k))
        target_error = abs(float(q[-1].sum()-min(mu.sum()/2, caps.sum())))
        assert error < 1e-9 and target_error < 1e-9
        results.append({"case": name, "means": means, "caps": limits,
                        "max_root_error_MW": error, "retuned_target_error_MW": target_error})
    return results


def main():
    source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    study = load_source()
    assert np.array_equal(study.W, PROBABILITIES) and np.array_equal(study.Z, STATES)
    cases = [audit_state(study, study.SEED+731, fractions, step, scale)
             for fractions in [(0.1, 0.1), (0.5, 0.5), (0.9, 0.9), (0.1, 0.9)]
             for step in [0, 24, 48, 76, 95]
             for scale in [0., .5, 1., 2.]]
    edges = audit_solver_edges(study)
    assert source_hash == hashlib.sha256(SOURCE.read_bytes()).hexdigest(), "source changed during audit"
    metrics = [k for k in cases[0] if k.startswith("max_") or k == "forecast_load_error"]
    report = {
        "status": "passed", "source_sha256": source_hash,
        "audit_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": {
            "n_operators": 2, "conditional_states_audited": len(cases),
            "full_joint_scenarios_per_state": 81,
            "total_joint_realizations_enumerated": 81*len(cases),
            "independent_checks": [
                "direct joint settlement expectation including imports and throughput cost",
                "marginal export means and quantile caps from complete joint enumeration",
                "unilateral quadratic objective fitted from direct settlement and price-spread evaluations",
                "entire controller objective multiplied by interval duration",
                "independent bisection versus sorted-breakpoint equilibrium solver",
                "conditional shortfall probability, power balance, and SOC recurrence",
            ],
            "limits": "Numerical finite-case audit, not a general proof or a repeat of trajectory experiments. Does not assess communication or transient electrical dynamics.",
        },
        "maxima": {key: max(c[key] for c in cases) for key in metrics},
        "max_import_in_enumerated_states_MW": max(c["maximum_import_MW"] for c in cases),
        "solver_edge_cases": edges, "conditional_cases": cases,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"status": report["status"], "cases": len(cases),
                      "maxima": report["maxima"], "output": str(OUTPUT)}, indent=2))


if __name__ == "__main__":
    main()
