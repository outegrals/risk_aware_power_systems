#!/usr/bin/env python3
"""Reproducible analytical/numerical checks for uncertainty and capped bidding.

Standard library only. This script does not alter the manuscript or MATLAB.
Forecasts x, bid quantities q, caps u, aggregate X=sum(x), Q=sum(q), and
coefficients k_i=1+4*a2*beta_i >= 0 use the manuscript's single-period model.

Uncertainty checks use direct settlement evaluation at an 81-node product
Gaussian quadrature rule, exact for the degree-four polynomials involved.
The analytical Gaussian variance is for realized revenue, not the squared
day-ahead/real-time price spread. No physical dispatch simulation is claimed.

The capped equilibrium solves r + sum(min(u_i,k_i*r)) = X. A separately
evaluated unilateral quadratic best response and KKT residual validate it.
Coordinated coefficient retuning is distinct from equilibrium at fixed k.
"""

from __future__ import annotations

import argparse
from itertools import product
import json
import math
from pathlib import Path
import random


SEED = 20261118
A2 = 0.1
A1 = 1.0
LOAD = 400.0


def capped_equilibrium(x: list[float], u: list[float], k: list[float]) -> dict:
    assert len(x) == len(u) == len(k)
    assert all(value > 0 for value in x)
    assert all(value >= 0 and math.isfinite(value) for value in u + k)
    total = sum(x)
    lower, upper = 0.0, total
    for _ in range(100):
        residual = (lower + upper) / 2
        if residual + sum(min(cap, coeff*residual) for cap, coeff in zip(u, k)) > total:
            upper = residual
        else:
            lower = residual
    residual = (lower + upper) / 2
    q = [min(cap, coeff*residual) for cap, coeff in zip(u, k)]
    return {"q": q, "Q": sum(q), "r": residual,
            "scalar_equation_residual": residual + sum(q) - total}


def unconstrained_equilibrium(x: list[float], k: list[float]) -> list[float]:
    residual = sum(x) / (1 + sum(k))
    return [coeff * residual for coeff in k]


def retuned_common_coefficient(x: list[float], u: list[float]) -> dict:
    total = sum(x)
    target = min(total/2, sum(u))
    if sum(u) <= total/2:
        level = max(u, default=0.0)
    else:
        lower, upper = 0.0, max(u)
        for _ in range(100):
            level = (lower + upper) / 2
            if sum(min(cap, level) for cap in u) > target:
                upper = level
            else:
                lower = level
        level = (lower + upper) / 2
    coeff = level / (total-target)
    return {"target_Q": target, "water_level_t": level, "common_k": coeff,
            "common_beta": (coeff-1)/(4*A2)}


def bid_objective(q_i: float, other_total: float, total: float, k_i: float) -> float:
    """J_i with bid-independent stochastic and revenue constants removed."""
    remaining = total-other_total-q_i
    return 2*A2*q_i*remaining - A2*(k_i-1)*remaining**2


def validate_equilibrium(x: list[float], u: list[float], k: list[float],
                         q: list[float], require_equilibrium: bool = True) -> dict:
    total, bid_total = sum(x), sum(q)
    br, gains, stable_gains, kkt = [], [], [], []
    for index, (cap, coeff, bid) in enumerate(zip(u, k, q)):
        other_total = sum(q[:index]) + sum(q[index+1:])
        response = max(0.0, min(cap, coeff/(1+coeff)*(total-other_total)))
        br.append(response)
        gains.append(bid_objective(response, other_total, total, coeff)
                     - bid_objective(bid, other_total, total, coeff))
        derivative_bracket = coeff*(total-bid_total)-bid
        displacement = response-bid
        stable_gains.append(2*A2*(derivative_bracket*displacement
                                 - 0.5*(1+coeff)*displacement**2))
        tolerance = 1e-11*max(1.0, total)
        if cap == 0:
            kkt.append(0.0)
        elif bid <= tolerance:
            kkt.append(max(0.0, derivative_bracket))
        elif cap-bid <= tolerance:
            kkt.append(max(0.0, -derivative_bracket))
        else:
            kkt.append(abs(derivative_bracket))
    violation = max([0.0] + [-bid for bid in q] + [bid-cap for bid, cap in zip(q, u)])
    metrics = {
        "direct_best_responses": br,
        "max_bid_feasibility_violation": violation,
        "max_best_response_distance": max(abs(a-b) for a, b in zip(br, q)),
        "max_kkt_bracket_residual": max(kkt),
        "max_normalized_kkt_residual": max(kkt)/max(1.0, total),
        "max_direct_unilateral_objective_gain": max([0.0]+gains),
        "max_stable_unilateral_objective_gain": max([0.0]+stable_gains),
    }
    if require_equilibrium:
        assert violation < 1e-10*max(1.0, total)
        assert metrics["max_best_response_distance"] < 1e-10*max(1.0, total)
        assert metrics["max_normalized_kkt_residual"] < 1e-10
        assert metrics["max_direct_unilateral_objective_gain"] < 1e-10*max(1.0, A2*total**2)
    return metrics


def deterministic_profit(x: list[float], q: list[float]) -> list[float]:
    total, bid_total = sum(x), sum(q)
    price_mean = A1+2*A2*(LOAD-total)
    return [price_mean*forecast + 2*A2*bid*(total-bid_total)
            for forecast, bid in zip(x, q)]


def known_cap_cases() -> list[dict]:
    configurations = [
        ("S1", [120.0, 100.0, 80.0], None),
        ("S2", [120.0, 100.0, 30.0], None),
        ("S3", [210.0, 45.0, 45.0], None),
        ("midpoint_counterexample", [1.0, 18.0, 18.0, 18.0, 45.0], None),
        ("shortage_caps", [120.0, 100.0, 80.0], [30.0, 20.0, 10.0]),
        ("zero_caps", [120.0, 100.0, 80.0], [0.0, 0.0, 0.0]),
    ]
    cases = []
    for name, x, selected_caps in configurations:
        u = x.copy() if selected_caps is None else selected_caps
        design = retuned_common_coefficient(x, u)
        n, total = len(x), sum(x)
        targets = [("risk_neutral_fixed", [1.0]*n),
                   ("ers_original_coefficient_fixed", [1/n]*n),
                   ("coordinated_retuned_common", [design["common_k"]]*n)]
        strategies = {}
        for label, k in targets:
            solution = capped_equilibrium(x, u, k)
            q = solution["q"]
            profits = deterministic_profit(x, q)
            optimum = sum((A1+2*A2*(LOAD-total))*v for v in x)
            optimum += 2*A2*design["target_Q"]*(total-design["target_Q"])
            strategies[label] = {
                **solution, "k": k, "beta": [(value-1)/(4*A2) for value in k],
                "alpha": [bid/forecast for bid, forecast in zip(q, x)],
                "individual_expected_profit": profits,
                "aggregate_expected_profit": sum(profits),
                "aggregate_profit_gap_to_capped_optimum": optimum-sum(profits),
                "validation": validate_equilibrium(x, u, k, q),
            }
        tuned = strategies["coordinated_retuned_common"]
        assert abs(tuned["Q"]-design["target_Q"]) < 1e-10*total
        assert abs(tuned["aggregate_profit_gap_to_capped_optimum"]) < 1e-9*max(1.0,total)
        # Deliberately evaluate the old componentwise clipping shortcut.
        clipped = {}
        for label, k in targets[:2]:
            q = [min(cap, bid) for cap, bid in zip(u, unconstrained_equilibrium(x, k))]
            clipped[label] = {"q": q, "Q": sum(q),
                              "validation": validate_equilibrium(x, u, k, q, False)}
        case = {"name": name, "x": x, "u": u, "X": total, "U": sum(u),
                "profit_assumption": "zero covariance; deterministic expected-profit table",
                "retuning_design": design, "strategies": strategies,
                "naive_componentwise_clipping": clipped}
        if name == "midpoint_counterexample":
            cutoff = (min(x)+max(x))/2
            saturated = [i for i, value in enumerate(x) if value < cutoff]
            remaining = n-len(saturated)
            level = (total/2-sum(x[i] for i in saturated))/remaining
            old_q = [value if i in saturated else level for i, value in enumerate(x)]
            assert old_q == [1.0, 18.0, 18.0, 18.0, -5.0]
            case["old_midpoint_rule"] = {"cutoff_in_mw": cutoff,
                "saturated_indices_zero_based": saturated,
                "saturated_total": sum(x[i] for i in saturated),
                "q": old_q, "alpha": [bid/forecast for bid,forecast in zip(old_q,x)],
                "contains_negative_bid": True}
        cases.append(case)
    by_name = {case["name"]: case for case in cases}
    expected_q = {
        "S1": {"risk_neutral_fixed": [75,75,75],
               "ers_original_coefficient_fixed": [50,50,50],
               "coordinated_retuned_common": [50,50,50]},
        "S2": {"risk_neutral_fixed": [220/3,220/3,30],
               "ers_original_coefficient_fixed": [44,44,30],
               "coordinated_retuned_common": [47.5,47.5,30]},
        "S3": {"risk_neutral_fixed": [105,45,45],
               "ers_original_coefficient_fixed": [52.5,45,45],
               "coordinated_retuned_common": [60,45,45]},
        "midpoint_counterexample": {"coordinated_retuned_common": [1,12.25,12.25,12.25,12.25]},
        "shortage_caps": {"coordinated_retuned_common": [30,20,10]},
    }
    for name, strategies in expected_q.items():
        for strategy, expected in strategies.items():
            actual = by_name[name]["strategies"][strategy]["q"]
            assert max(abs(a-b) for a,b in zip(actual,expected)) < 1e-9
    assert by_name["S2"]["naive_componentwise_clipping"]["risk_neutral_fixed"]["validation"]["max_direct_unilateral_objective_gain"] > 0
    return cases


def cholesky(matrix: list[list[float]]) -> list[list[float]]:
    n = len(matrix)
    lower = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            value = matrix[i][j]-sum(lower[i][h]*lower[j][h] for h in range(j))
            if i == j:
                assert value > 0, "Base covariance is not positive definite"
                lower[i][j] = math.sqrt(value)
            else:
                lower[i][j] = value/lower[j][j]
    return lower


def base_covariance() -> list[list[float]]:
    std = [12.0, 10.0, 8.0, 20.0]
    return [[std[i]*std[j]*(1.0 if i == j else (0.2 if 3 in (i,j) else 0.3))
             for j in range(4)] for i in range(4)]


def gaussian_revenue_moments(x: list[float], q: list[float], covariance: list[list[float]]) -> dict:
    n, total, bid_total = len(x), sum(x), sum(q)
    variance_e = covariance[n][n]-2*sum(covariance[i][n] for i in range(n))
    variance_e += sum(covariance[i][j] for i in range(n) for j in range(n))
    variance_price = 4*A2**2*variance_e
    mean_price = A1+2*A2*(LOAD-total)
    expected, variance, corrections = [], [], []
    for i in range(n):
        cov_price_gen = 2*A2*(covariance[i][n]-sum(covariance[i][:n]))
        corrections.append(cov_price_gen)
        expected.append(mean_price*x[i]+cov_price_gen+2*A2*q[i]*(total-bid_total))
        exposure = x[i]-q[i]
        variance.append(mean_price**2*covariance[i][i]+exposure**2*variance_price
                        +2*mean_price*exposure*cov_price_gen
                        +variance_price*covariance[i][i]+cov_price_gen**2)
    return {"V_net_load_error": variance_e, "variance_real_time_price": variance_price,
            "mean_real_time_price": mean_price,
            "mean_day_ahead_price": A1+2*A2*(LOAD-bid_total),
            "price_spread_squared_bias": 4*A2**2*(total-bid_total)**2,
            "price_spread_second_moment_R": variance_price+4*A2**2*(total-bid_total)**2,
            "individual_covariance_profit_corrections": corrections,
            "aggregate_covariance_profit_correction": sum(corrections),
            "individual_expected_revenue": expected, "aggregate_expected_revenue": sum(expected),
            "individual_realized_revenue_variance": variance,
            "individual_realized_revenue_standard_deviation": [math.sqrt(max(0,v)) for v in variance]}


def direct_gaussian_quadrature(x: list[float], q: list[float], lower: list[list[float]], scale: float) -> dict:
    nodes = [(-math.sqrt(3), 1/6), (0.0, 2/3), (math.sqrt(3), 1/6)]
    total, bid_total = sum(x), sum(q)
    day_ahead_price = A1+2*A2*(LOAD-bid_total)
    observations = []
    for points in product(nodes, repeat=4):
        weight = math.prod(point[1] for point in points)
        errors = [scale*sum(lower[i][j]*points[j][0] for j in range(4)) for i in range(4)]
        generation = [value+error for value,error in zip(x,errors[:3])]
        load = LOAD+errors[3]
        real_time_price = A1+2*A2*(load-sum(generation))
        revenue = [day_ahead_price*bid+real_time_price*(actual-bid)
                   for actual,bid in zip(generation,q)]
        observations.append((weight,revenue,(real_time_price-day_ahead_price)**2))
    means = [sum(weight*revenue[i] for weight,revenue,_ in observations) for i in range(3)]
    variances = [sum(weight*(revenue[i]-means[i])**2 for weight,revenue,_ in observations) for i in range(3)]
    risk = sum(weight*spread for weight,_,spread in observations)
    return {"n_quadrature_nodes": len(observations), "mean_revenue": means,
            "variance_revenue": variances, "R": risk}


def uncertainty_sensitivity() -> dict:
    x = [120.0,100.0,80.0]
    base = base_covariance()
    lower = cholesky(base)
    reconstructed = [[sum(lower[i][k]*lower[j][k] for k in range(4)) for j in range(4)] for i in range(4)]
    reconstruction_error = max(abs(base[i][j]-reconstructed[i][j]) for i in range(4) for j in range(4))
    assert reconstruction_error < 1e-10
    coefficients = {"risk_neutral": [1.0]*3, "ers": [1/3]*3}
    profiles = {name:unconstrained_equilibrium(x,k) for name,k in coefficients.items()}
    records = []
    for scale in (0.0,0.5,1.0,2.0):
        covariance = [[scale**2*value for value in row] for row in base]
        for strategy,q in profiles.items():
            analytic = gaussian_revenue_moments(x,q,covariance)
            independent = direct_gaussian_quadrature(x,q,lower,scale)
            mean_error = max(abs(a-b) for a,b in zip(analytic["individual_expected_revenue"],independent["mean_revenue"]))
            variance_error = max(abs(a-b) for a,b in zip(analytic["individual_realized_revenue_variance"],independent["variance_revenue"]))
            risk_error = abs(analytic["price_spread_second_moment_R"]-independent["R"])
            # Independently fit each unilateral objective's quadratic from
            # direct settlement quadrature at three feasible bids. The fitted
            # maximizers verify that covariance scaling leaves bids unchanged.
            fitted_best_responses = []
            for i in range(3):
                step = x[i]/2
                beta = (coefficients[strategy][i]-1)/(4*A2)
                values = []
                for candidate in (0.0,step,2*step):
                    trial_q = q.copy()
                    trial_q[i] = candidate
                    moments = direct_gaussian_quadrature(x,trial_q,lower,scale)
                    values.append(moments["mean_revenue"][i]-beta*moments["R"])
                quadratic = (values[2]-2*values[1]+values[0])/(2*step**2)
                linear = (values[1]-values[0]-quadratic*step**2)/step
                assert quadratic < 0
                fitted_best_responses.append(max(0.0,min(x[i],-linear/(2*quadratic))))
            best_response_error = max(abs(a-b) for a,b in zip(q,fitted_best_responses))
            assert mean_error < 1e-8
            assert variance_error < 1e-7
            assert risk_error < 1e-9
            assert best_response_error < 1e-8
            records.append({"uncertainty_std_scale":scale, "covariance_multiplier":scale**2,
                "strategy":strategy,"q":q,"alpha":[bid/forecast for bid,forecast in zip(q,x)],
                **analytic,"direct_settlement_quadrature":independent,
                "quadrature_fitted_unilateral_best_responses":fitted_best_responses,
                "validation":{"max_mean_error":mean_error,"max_variance_error":variance_error,
                              "R_error":risk_error,"max_quadrature_best_response_error":best_response_error}})
    return {"x":x,"load_forecast":LOAD,"a2":A2,"a1":A1,"base_covariance":base,
            "base_covariance_positive_definite":True,"cholesky_diagonal":[lower[i][i] for i in range(4)],
            "cholesky_reconstruction_max_error":reconstruction_error,
            "method":"81-node product Gaussian quadrature, exact through degree five in each independent normal coordinate",
            "distribution_scope":"untruncated joint Gaussian; no resampling or clipping; quantities are market settlement revenues",
            "variance_identity":"Var(revenue_i)=m_lambda^2*sigma_ii+(x_i-q_i)^2*Var(lambda)+2*m_lambda*(x_i-q_i)*Cov(lambda,x_i)+Var(lambda)*sigma_ii+Cov(lambda,x_i)^2",
            "records":records}


def random_cap_checks() -> dict:
    rng = random.Random(SEED)
    records = []
    for n in (2,3,10,50,100):
        for replicate in range(20):
            x = [rng.uniform(1,200) for _ in range(n)]
            # Includes zero caps, aggregate shortages, and ordinary forecast caps.
            mode = replicate % 4
            u = [(0.0 if (mode == 0 and i % 4 == 0) else
                  value*(rng.uniform(0.01,0.2) if mode == 1 else rng.uniform(0.2,1.0)))
                 for i,value in enumerate(x)]
            k = [(0.0 if i % 4 == 0 else 10**rng.uniform(-3,1)) for i in range(n)]
            # A distinct all-zero-coefficient instance exercises beta=-1/(4*a2).
            if replicate == 0:
                k = [0.0]*n
            fixed = capped_equilibrium(x,u,k)
            fixed_validation = validate_equilibrium(x,u,k,fixed["q"])
            design = retuned_common_coefficient(x,u)
            tuned_k = [design["common_k"]]*n
            tuned = capped_equilibrium(x,u,tuned_k)
            tuned_validation = validate_equilibrium(x,u,tuned_k,tuned["q"])
            target_error = abs(tuned["Q"]-design["target_Q"])
            assert target_error < 1e-10*sum(x)
            records.append({"n":n,"replicate":replicate,"x":x,"u":u,"k_fixed":k,
                "fixed_equilibrium":fixed,"fixed_validation":fixed_validation,
                "retuning_design":design,"retuned_equilibrium":tuned,
                "retuned_validation":tuned_validation,"target_Q_absolute_error":target_error})
    validations = [record[key] for record in records for key in ("fixed_validation","retuned_validation")]
    return {"seed":SEED,"number_of_instances":len(records),"equilibria_checked":len(validations),
            "n_values":[2,3,10,50,100],"replicates_per_n":20,
            "max_normalized_kkt_residual":max(v["max_normalized_kkt_residual"] for v in validations),
            "max_best_response_distance":max(v["max_best_response_distance"] for v in validations),
            "max_direct_unilateral_objective_gain":max(v["max_direct_unilateral_objective_gain"] for v in validations),
            "max_stable_unilateral_objective_gain":max(v["max_stable_unilateral_objective_gain"] for v in validations),
            "max_retuned_target_Q_absolute_error":max(record["target_Q_absolute_error"] for record in records),
            "records":records}


def vector_text(values: list[float]) -> str:
    return "("+", ".join(f"{value:.6g}" for value in values)+")"


def markdown_report(risk: dict, caps: list[dict], random_checks: dict) -> str:
    lines = ["# Risk and capped-equilibrium validation", "",
        "Generated by `python3 scripts/validate_risk_and_caps.py`; no simulations of grid dynamics or strategic beta selection are implied.", "",
        "## Uncertainty sensitivity", "",
        "S1 uses x=(120,100,80) MW, load=400 MW, a2=0.1, a1=1. Covariance is scaled by t² while the means are fixed. All profit entries include the covariance correction. Price variance and the squared mean price spread are reported separately.", "",
        "| t | Strategy | Alpha | Var(real-time price) | Squared price-spread bias | R | Total expected revenue | Individual revenue SD |",
        "|---:|---|---|---:|---:|---:|---:|---|"]
    for row in risk["records"]:
        lines.append(f"| {row['uncertainty_std_scale']:g} | {row['strategy']} | {vector_text(row['alpha'])} | {row['variance_real_time_price']:.6f} | {row['price_spread_squared_bias']:.6f} | {row['price_spread_second_moment_R']:.6f} | {row['aggregate_expected_revenue']:.6f} | {vector_text(row['individual_realized_revenue_standard_deviation'])} |")
    lines += ["", "The invariant bids follow from the stated model. ERS changes the squared mean price spread while real-time price variance is identical across strategies. An 81-node product Gaussian quadrature independently evaluates settlement revenues and their variances; it agrees with the analytical expressions within the tolerances recorded in JSON. Unilateral quadratics fitted from three direct quadrature evaluations also recover the same best-response bids at every uncertainty scale. ERS increases every operator's realized-revenue standard deviation in this example, so its expected-profit improvement does not demonstrate variance reduction.", "",
        "## Bounded-game equilibria", "",
        "This table uses zero covariance to isolate the capacity correction and reproduce the original deterministic profit convention. Its profits must not be mixed with the stochastic totals above. Each fixed row solves the bounded game at its stated coefficient; the retuned row explicitly changes the common coefficient by coordination.", "",
        "| Case | Strategy | Common k | Common beta | q (MW) | Q | Total expected profit | Gap to capped aggregate optimum |",
        "|---|---|---:|---:|---|---:|---:|---:|"]
    for case in caps:
        for label,row in case["strategies"].items():
            lines.append(f"| {case['name']} | {label} | {row['k'][0]:.6g} | {row['beta'][0]:.6g} | {vector_text(row['q'])} | {row['Q']:.6f} | {row['aggregate_expected_profit']:.6f} | {row['aggregate_profit_gap_to_capped_optimum']:.6f} |")
    lines += ["", "The old midpoint rule on x=(1,18,18,18,45) assigns q=(1,18,18,18,-5), so it is infeasible. The retuned construction assigns q=(1,12.25,12.25,12.25,12.25). With shortage caps u=(30,20,10), the aggregate optimum is Q=60 rather than X/2=150; the construction supports all caps with common k=1/8.", "",
        "## Independent best-response and randomized checks", "",
        f"Seed {SEED}; {random_checks['number_of_instances']} capacity/coefficient instances at n=2,3,10,50,100; {random_checks['equilibria_checked']} fixed or retuned equilibria. Cases include k=0, zero caps, aggregate shortages, and heterogeneous nonzero coefficients.", "",
        f"- Maximum normalized KKT bracket residual: {random_checks['max_normalized_kkt_residual']:.6g}.",
        f"- Maximum distance to direct unilateral best response: {random_checks['max_best_response_distance']:.6g} MW.",
        f"- Maximum direct objective gain at the computed best response: {random_checks['max_direct_unilateral_objective_gain']:.6g} (floating-point cancellation scale).",
        f"- Maximum gain using the stable quadratic-difference expression: {random_checks['max_stable_unilateral_objective_gain']:.6g}.",
        f"- Maximum retuned aggregate-target error: {random_checks['max_retuned_target_Q_absolute_error']:.6g} MW.", "",
        "These are reproducible finite numerical checks alongside the formulas, not proofs of the general theorem, performance benchmarks, or evidence for dynamic microgrids. They do not establish incentive compatibility when beta itself is a strategic choice. General-cap computation may disclose quantities outside the separate masked-uncapped privacy theorem.", "",
        "No plotting package is required or invoked. The JSON contains all table values for plotting.", ""]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,default=Path(__file__).resolve().parents[1]/"output/risk_caps")
    args = parser.parse_args()
    risk = uncertainty_sensitivity()
    caps = known_cap_cases()
    random_checks = random_cap_checks()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    for name,value in [("uncertainty_sensitivity.json",risk),("known_cap_cases.json",caps),
                       ("random_cap_checks.json",random_checks)]:
        (args.output_dir/name).write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")
    (args.output_dir/"validation_report.md").write_text(markdown_report(risk,caps,random_checks))
    print(f"All checks passed; wrote evidence to {args.output_dir}")
    baseline = next(row for row in risk["records"] if row["uncertainty_std_scale"] == 1)
    print(f"Baseline covariance V={baseline['V_net_load_error']}; {random_checks['equilibria_checked']} equilibria checked")
    print(f"Maximum normalized KKT residual: {random_checks['max_normalized_kkt_residual']:.6g}")
    print(f"Maximum stable best-response gain: {random_checks['max_stable_unilateral_objective_gain']:.6g}")


if __name__ == "__main__":
    main()
