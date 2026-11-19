#!/usr/bin/env python3
"""Reproducible masked-consensus communication and convergence study.

This is a numerical/network model, not a distributed deployment or a
cryptographic implementation. Public seeded masks must never be production
secrets. Run with Python 3 and NumPy; output defaults to output/communication.
"""

from __future__ import annotations

import argparse
from collections import Counter
from decimal import Decimal, localcontext, ROUND_CEILING
from fractions import Fraction
import json
import math
from pathlib import Path
import platform
import random
import sys
import time

import numpy as np

from validate_masked_consensus import masked_values, next_prime

BASE_SEED = 20261118
GRAPHS = ("star", "ring", "circulant4", "complete")


def graph(n: int, kind: str) -> dict:
    m = n + 1
    if kind == "star":
        edges = [(i, n) for i in range(n)]
        shifts = []
    elif kind == "complete":
        edges = [(i, j) for i in range(m) for j in range(i + 1, m)]
        shifts = []
    else:
        step = max(2, math.isqrt(m))
        shifts = sorted({1, m - 1} if kind == "ring" else
                        {1, m - 1, step % m, (-step) % m} - {0})
        edges = sorted({tuple(sorted((i, (i + s) % m)))
                        for i in range(m) for s in shifts})
    degrees = [0] * m
    for i, j in edges:
        degrees[i] += 1
        degrees[j] += 1
    dmax = max(degrees)
    denominator = dmax + 1
    if kind == "star":
        rho_upper = Fraction(n, m)
        rho = float(rho_upper)
    elif kind == "complete":
        rho_upper = Fraction(0)
        rho = 0.0
    else:
        # The exact circulant Laplacian eigenvalues are
        # 2 sum_s sin^2(pi*j*s/m). For distance t to the nearest integer,
        # 0 <= t <= 1/2 and sin(pi*t) >= 2*t. Thus lambda_j is at least
        # 8 sum_s dist(j*s,m Z)^2 / m^2. All arithmetic below is exact.
        lower_gap = Fraction(8 * min(
            sum(min((j*s) % m, m - (j*s) % m)**2 for s in shifts)
            for j in range(1, m)), m*m)
        rho_upper = max(1 - lower_gap / denominator,
                        Fraction(2*dmax, denominator) - 1)
        eigenvalues = [1 - sum(1 - math.cos(2*math.pi*j*s/m)
                               for s in shifts) / denominator
                       for j in range(1, m)]
        rho = max(abs(x) for x in eigenvalues)
        assert rho <= float(rho_upper) + 1e-12
    # Greedy edge coloring produces matchings. During a color's slot, each
    # node exchanges a pair of directed scalar messages with at most one peer.
    used = [set() for _ in range(m)]
    color_counts: Counter = Counter()
    for i, j in edges:
        color = 0
        while color in used[i] or color in used[j]:
            color += 1
        used[i].add(color)
        used[j].add(color)
        color_counts[color] += 1
    return {"kind": kind, "m": m, "edges": edges, "shifts": shifts,
            "degrees": degrees, "dmax": dmax, "denominator": denominator,
            "rho_numerical": rho, "rho_upper": rho_upper,
            "color_directed_messages": [2*color_counts[c]
                                        for c in sorted(color_counts)]}


def step(z: np.ndarray, spec: dict) -> np.ndarray:
    """Positive weighted sum implementation of W=I-L/(dmax+1)."""
    n, kind = spec["m"] - 1, spec["kind"]
    if kind == "star":
        ans = np.empty_like(z)
        ans[:-1] = (n*z[:-1] + z[-1]) / (n+1)
        ans[-1] = (np.sum(z[:-1], axis=0) + z[-1]) / (n+1)
        return ans
    if kind == "complete":
        return np.broadcast_to(np.sum(z, axis=0) / (n+1), z.shape).copy()
    ans = z.copy()
    for shift in spec["shifts"]:
        ans += np.roll(z, shift, axis=0)
    return ans / spec["denominator"]


def upward_power(x: Decimal, k: int) -> Decimal:
    """Multiplication-only positive integer power, rounded upward by context."""
    result = Decimal(1)
    while k:
        if k & 1:
            result *= x
        k //= 2
        if k:
            x *= x
    return result


def certified_budget(spec: dict, prime: int) -> dict:
    m, rho = spec["m"], spec["rho_upper"]
    root_ceiling = math.isqrt(m)
    root_ceiling += root_ceiling * root_ceiling < m
    constant = m * root_ceiling * (prime - 1)
    k = 1 if not rho else max(1, math.ceil(
        math.log(0.20 / constant) / math.log(float(rho))))
    with localcontext() as ctx:
        ctx.prec = 70
        ctx.rounding = ROUND_CEILING
        r = Decimal(rho.numerator) / Decimal(rho.denominator)
        ideal = Decimal(constant) * upward_power(r, k)
        while ideal >= Decimal("0.25"):
            k += 1
            ideal = Decimal(constant) * upward_power(r, k)
        # Every row uses positive integer-weighted sums and one division.
        # At most m+5 rounded operations suffice for any row in any kernel.
        # Under IEEE binary64, RN, no overflow/underflow, each local relative
        # error <= gamma_t. W is nonnegative with row sum one. Induction gives
        # ||z_hat(k)-W^k z(0)||_infty <= M*((1+gamma_t)^k-1).
        # Include final scaling by m, and use exact input bound M=p-1.
        u = Decimal(1) / Decimal(2**53)
        t = m + 5
        gamma = Decimal(t)*u / (1 - Decimal(t)*u)
        roundoff = Decimal(m*(prime - 1)) * (
            upward_power(1 + gamma, k) * (1 + u) - 1)
        total = ideal + roundoff
        assert total < Decimal("0.5"), (spec["kind"], m, total)
    return {"public_iterations": k,
            "ideal_scaled_error_bound_upward": str(ideal),
            "binary64_scaled_error_bound_upward": str(roundoff),
            "combined_scaled_error_bound_upward": str(total),
            "rounding_certified_under_stated_arithmetic_model": True,
            "row_operation_budget": t}


def inputs(n: int, graph_id: int, trials: int, bound: int = 120,
           delta: float = 1.0) -> dict:
    forecasts = ([120, 100, 80] if n == 3 and bound == 120 else
                 [random.Random(BASE_SEED+n*1000+i).randint(2*bound//3, bound)
                  for i in range(n)])
    prime = 1009 if n == 3 and bound == 120 else next_prime(n*bound+1)
    mask_edges = [(i, (i+1) % n) for i in range(n)]
    runs = []
    for trial in range(trials):
        seed = BASE_SEED + 100000*n + 1000*graph_id + trial + 10*bound
        rng = random.Random(seed)
        masks = [rng.randrange(prime) for _ in mask_edges]
        values = masked_values(forecasts, mask_edges, masks, prime)
        assert sum(values) % prime == sum(forecasts)
        runs.append({"seed": seed, "public_test_masks": masks,
                     "masked_values": values,
                     "residue_sum": sum(values),
                     "modular_wraps": (sum(values)-sum(forecasts))//prime})
    return {"forecasts_encoded": forecasts, "forecasts_MW": [delta*x for x in forecasts],
            "total": sum(forecasts), "prime": prime, "B": bound, "delta_MW": delta,
            "mask_edges": mask_edges, "runs": runs}


def spectral_state(initial: np.ndarray, spec: dict, k: int) -> np.ndarray:
    """Fast evaluation of ideal W^k; explicitly not a full iteration execution."""
    if k == 0:
        return initial.copy()
    mean = np.mean(initial, axis=0)
    if spec["kind"] == "complete":
        return np.broadcast_to(mean, initial.shape).copy()
    if spec["kind"] == "star":
        z = np.broadcast_to(mean, initial.shape).copy()
        z[:-1] += ((spec["m"]-1)/spec["m"])**k * (
            initial[:-1]-np.mean(initial[:-1], axis=0))
        return z
    m = spec["m"]
    eigen = 1 - sum(1-np.cos(2*np.pi*np.arange(m)*s/m)
                    for s in spec["shifts"])/spec["denominator"]
    return np.fft.ifft(np.fft.fft(initial, axis=0)*eigen[:, None]**k,
                       axis=0).real


def convergence(spec: dict, data: dict, budget: dict) -> dict:
    m, k = spec["m"], budget["public_iterations"]
    z = np.array([run["masked_values"] + [0] for run in data["runs"]],
                 dtype=np.float64).T
    initial = z.copy()
    totals = np.array([run["residue_sum"] for run in data["runs"]])
    first_correct = np.full(len(totals), -1, dtype=int)
    first_quarter = first_correct.copy()
    sample_rounds = set(np.unique(np.geomspace(1, k, min(90, k)).astype(int)))
    sample_rounds.update([0, k])
    trace = []
    max_sum_drift = 0.0
    started = time.perf_counter()
    full_execution = k <= 12000
    iteration_indices = range(k+1) if full_execution else sorted(sample_rounds)
    for iteration in iteration_indices:
        if not full_execution:
            z = spectral_state(initial, spec, iteration)
        errors = np.max(np.abs(m*z - totals), axis=0)
        if np.any(first_correct < 0):
            correct = np.all(np.rint(m*z).astype(np.int64) == totals, axis=0)
            first_correct[(first_correct < 0) & correct] = iteration
        first_quarter[(first_quarter < 0) & (errors < .25)] = iteration
        if iteration in sample_rounds:
            max_sum_drift = max(max_sum_drift,
                               float(np.max(np.abs(np.sum(z, axis=0)-totals))))
            trace.append({"iteration": int(iteration),
                          "max_scaled_error_over_trials": float(max(errors))})
        if full_execution and iteration < k:
            z = step(z, spec)
    if not full_execution:
        # Binary search the monotone max-norm error threshold; this uses an
        # offline oracle and never changes the public protocol budget.
        for trial in range(len(totals)):
            for threshold, destination in ((.5, first_correct), (.25, first_quarter)):
                lo, hi = 0, k
                while lo < hi:
                    mid = (lo+hi)//2
                    candidate = spectral_state(initial[:, trial:trial+1], spec, mid)
                    error = np.max(np.abs(m*candidate-totals[trial]))
                    if error < threshold:
                        hi = mid
                    else:
                        lo = mid+1
                destination[trial] = lo
    validation_seconds = time.perf_counter() - started
    final_errors = np.max(np.abs(m*z-totals), axis=0)
    decoded = np.rint(m*z).astype(np.int64) % data["prime"]
    assert np.all(decoded == data["total"])
    assert max(final_errors) < float(budget["combined_scaled_error_bound_upward"])
    assert np.all(first_correct >= 0) and np.all(first_quarter >= 0)
    return {"all_trials_all_nodes_decode_correctly": True,
            "evaluation_method": "full synchronous binary64 iteration" if full_execution else
                "spectral evaluation of ideal W^k; sampled states and offline threshold search",
            "final_max_scaled_errors": final_errors.tolist(),
            "offline_first_correct_decoding_rounds": first_correct.tolist(),
            "offline_first_below_quarter_rounds": first_quarter.tolist(),
            "sampled_max_total_drift_binary64": max_sum_drift,
            "validation_seconds_including_oracle_diagnostics": validation_seconds,
            "trace": trace, "initial": initial[:, 0]}


def kernel_timing(spec: dict, initial: np.ndarray, k: int, repeats: int) -> dict:
    # Warm up allocation/dispatch; time only unbatched consensus updates.
    z = initial.copy()
    for _ in range(min(k, 20)):
        z = step(z, spec)
    measured_updates = min(k, 2000)
    elapsed = []
    for _ in range(repeats):
        z = initial.copy()
        started = time.perf_counter()
        for _ in range(measured_updates):
            z = step(z, spec)
        elapsed.append(time.perf_counter()-started)
    spectral = spectral_state(initial[:, None], spec, measured_updates)[:, 0]
    comparison = float(np.max(np.abs(z-spectral)))
    # This comparison checks the fast diagnostic against executed updates;
    # it is separate from the conservative forward-error certificate.
    assert comparison <= 1e-7*max(1.0, float(np.max(initial)))
    return {"repeats": repeats, "seconds": elapsed,
            "measured_updates_per_repeat": measured_updates,
            "public_budget_executed_in_timing": measured_updates == k,
            "max_state_difference_from_spectral_at_measured_updates": comparison,
            "median_seconds": float(np.median(elapsed)),
            "minimum_seconds": min(elapsed), "maximum_seconds": max(elapsed),
            "excludes": ["mask exchange", "encryption", "network I/O",
                         "decoding", "diagnostic checks"],
            "interpretation": "one-host NumPy emulator, not distributed node latency"}


def max_geometric_pmf(messages: int, loss: float, batches: int) -> tuple:
    """Maximum retries in a batch; truncate with a documented union bound."""
    if loss == 0:
        return np.array([1.0]), 0.0, 1.0
    # P(max>A) <= messages*loss**A. Force total-run omitted mass < 1e-13.
    amax = max(1, math.ceil(math.log(1e-13/(messages*batches)) / math.log(loss)))
    cdf = np.array([math.exp(messages*math.log1p(-loss**a))
                    for a in range(1, amax+1)])
    pmf = np.diff(np.concatenate(([0.0], cdf)))
    omitted = messages * loss**amax * batches
    pmf[-1] += 1 - sum(pmf)
    expected = 1.0 + sum(-math.expm1(messages*math.log1p(-loss**a))
                         for a in range(1, amax+10))
    return pmf, omitted, expected


def network_model(spec: dict, k: int, loss: float, schedule: str,
                  episodes: int, seed: int, slot_ms: float = 1.0) -> dict:
    """Reliable buffered message delivery followed by one synchronous update.

    One slot is a hypothetical 1 ms data+ack attempt. Failed directed
    messages retry; successful ones do not. Within a batch all links run
    in parallel. Either every directed edge is one batch, or a greedy
    matching schedule restricts each node to one neighboring peer at once.
    States are frozen/buffered until all batches for the round finish.
    """
    batches = ([2*len(spec["edges"])] if schedule == "parallel_links" else
               spec["color_directed_messages"])
    n = spec["m"]-1
    if schedule == "parallel_links":
        setup_batches = [n]
    else:
        # A feasible coloring of the private operator cycle, one directed
        # field element per edge, separate from consensus scalar messages.
        setup_batches = ([n//2, n//2] if n % 2 == 0 else
                         [(n-1)//2, (n-1)//2, 1])
    group_counts = Counter()
    for count in batches:
        group_counts[count] += k
    for count in setup_batches:
        group_counts[count] += 1
    rng = np.random.default_rng(seed)
    slots = np.zeros(episodes, dtype=np.int64)
    expected_slots = 0.0
    truncation_probability_bound = 0.0
    for message_count, count in group_counts.items():
        pmf, omitted, expectation = max_geometric_pmf(message_count, loss, count)
        draws = rng.multinomial(count, pmf, size=episodes)
        slots += draws @ np.arange(1, len(pmf)+1)
        expected_slots += count*expectation
        truncation_probability_bound += omitted
    seconds = slots*slot_ms/1000
    logical = 2*len(spec["edges"])*k+n
    proportions = {str(d): float(np.mean(seconds <= d)) for d in (1, 10, 60, 300, 900)}
    intervals = {}
    z95 = 1.959963984540054
    for deadline, proportion in proportions.items():
        divisor = 1 + z95*z95/episodes
        center = (proportion+z95*z95/(2*episodes))/divisor
        radius = z95*math.sqrt(proportion*(1-proportion)/episodes+
                              z95*z95/(4*episodes*episodes))/divisor
        intervals[deadline] = [max(0, center-radius), min(1, center+radius)]
    return {"schedule": schedule, "independent_attempt_loss_probability": loss,
            "hypothetical_attempt_slot_ms": slot_ms, "seed": seed,
            "monte_carlo_episodes": episodes, "batches_per_round": len(batches),
            "mask_setup_batches": len(setup_batches),
            "expected_network_seconds": expected_slots*slot_ms/1000,
            "sample_mean_seconds": float(np.mean(seconds)),
            "sample_p50_seconds": float(np.quantile(seconds, .5)),
            "sample_p95_seconds": float(np.quantile(seconds, .95)),
            "sample_max_seconds": float(np.max(seconds)),
            "deadline_met_fractions": proportions,
            "deadline_fraction_wilson_95pct_intervals": intervals,
            "successful_directed_scalar_or_field_deliveries": logical,
            "expected_attempted_directed_scalar_transmissions": logical/(1-loss),
            "per_episode_tail_truncation_probability_upper_bound":
                truncation_probability_bound,
            "averaging_updates": k,
            "state_trace": "identical to lossless synchronous trace; only barriers delayed"}


def network_suite(spec: dict, k: int, episodes: int) -> list:
    n, graph_id = spec["m"]-1, GRAPHS.index(spec["kind"])
    return [network_model(spec, k, loss, schedule, episodes,
                          BASE_SEED+n*10000+graph_id*100+loss_id*10+schedule_id,
                          slot_ms=slot)
            for loss_id, loss in enumerate((0, .05, .2))
            for schedule_id, schedule in enumerate(("parallel_links", "one_peer_per_node"))
            for slot in (1, 10)]


def finite_retry_schedule_check() -> dict:
    """Event-level loss test checks the buffered implementation, not just times."""
    n, spec = 10, graph(10, "circulant4")
    data = inputs(n, 44, 1)
    z = np.array(data["runs"][0]["masked_values"]+[0], dtype=float)
    reference = z.copy()
    rng = np.random.default_rng(BASE_SEED+55)
    failed, attempts, slots = 0, 0, 0
    for _ in range(120):
        pending = {(i, j) for edge in spec["edges"] for i, j in (edge, edge[::-1])}
        received = {}
        while pending:
            slots += 1
            for i, j in sorted(pending):
                attempts += 1
                if rng.random() < .2:
                    failed += 1
                else:
                    received[(i, j)] = z[i]
                    pending.remove((i, j))
        updated = np.empty_like(z)
        for j in range(spec["m"]):
            # Same neighbor order as the array kernel for bitwise comparison.
            value = z[j]
            for shift in spec["shifts"]:
                value += received[((j-shift) % spec["m"], j)]
            updated[j] = value/spec["denominator"]
        reference = step(reference, spec)
        assert np.array_equal(updated, reference)
        z = updated
    return {"n_operators": n, "graph": spec["kind"], "rounds": 120,
            "loss_probability": .2, "attempted_scalar_messages": attempts,
            "failed_attempts": failed, "network_slots": slots,
            "every_round_bitwise_identical_to_lossless_reference": True,
            "seed": BASE_SEED+55}


def write_svg(cases: list, path: Path) -> None:
    chosen = [c for c in cases if c["n_operators"] == 100 and c["input_data"]["B"] == 10000]
    if not chosen:
        return
    colors = {"star": "#0969da", "ring": "#d1242f", "circulant4": "#1a7f37",
              "complete": "#8250df"}
    left, top, width, height = 92, 52, 680, 355
    xmax = math.log10(max(c["budget"]["public_iterations"] for c in chosen)+1)
    ymax = math.ceil(math.log10(max(t["max_scaled_error_over_trials"]
                               for c in chosen for t in c["convergence"]["trace"])))
    ymin = -3
    def xy(k, e):
        return (left+width*math.log10(k+1)/xmax,
                top+height*(ymax-math.log10(max(10**ymin, e)))/(ymax-ymin))
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="900" height="485" viewBox="0 0 900 485">',
           '<rect width="900" height="485" fill="white"/>',
           '<g font-family="Arial,sans-serif" font-size="13" fill="#222">',
           '<text x="92" y="25" font-size="18">Masked consensus: 100 operators, three public mask trials</text>']
    for exponent in range(ymin, ymax+1):
        _, y = xy(0, 10**exponent)
        svg += [f'<path d="M {left} {y:.2f} H {left+width}" stroke="#ddd"/>',
                f'<text x="{left-12}" y="{y+4:.2f}" text-anchor="end">10^{exponent}</text>']
    for exponent in range(math.floor(xmax)+1):
        x, _ = xy(10**exponent-1, 1)
        svg += [f'<text x="{x:.2f}" y="430" text-anchor="middle">{10**exponent-1}</text>']
    for index, case in enumerate(chosen):
        points = " ".join(f'{x:.2f},{y:.2f}' for x,y in
                          (xy(t["iteration"], t["max_scaled_error_over_trials"])
                           for t in case["convergence"]["trace"]))
        color = colors[case["graph"]]
        svg += [f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2"/>',
                f'<text x="{left+index*180}" y="458" fill="{color}">{case["graph"]}</text>']
    _, y = xy(0, .5)
    svg += [f'<path d="M {left} {y:.2f} H {left+width}" stroke="#555" stroke-dasharray="5,4"/>',
            f'<text x="{left+width+6}" y="{y+4:.2f}">0.5</text>',
            '<text x="432" y="480" text-anchor="middle">Consensus updates (log scale of k+1)</text>',
            '<text transform="translate(19,245) rotate(-90)" text-anchor="middle">Maximum scaled residue-sum error</text>',
            '</g></svg>']
    path.write_text("\n".join(svg)+"\n")


def write_report(report: dict, directory: Path) -> None:
    lines = ["# Communication and convergence validation", "",
             "Run `python3 scripts/validate_communication.py`. Python 3 and NumPy are sufficient.", "",
             "These are mathematical and numerical checks plus a synthetic network model. "
             "The timing measurements are from one host running all nodes. They are not deployed microgrid latency, "
             "secure-channel measurements, or evidence of physical closed-loop stability.", "",
             "## Protocol and guarantees", "",
             "- The public consensus graph contains n operators and the utility. The private mask graph is an operator cycle; its honest-induced connectivity still limits collusion tolerance. Its private logical links are assumed directly available with the same attempt-slot model; additional relay/routing costs are excluded.",
             "- The fine-grid profile uses public B=10000 integer units and Delta=.001 MW (forecasts at most 10 MW), at n=10,50,100. Its forecasts are independently generated; this is not a communication replay of the separate physical trajectory. A baseline profile uses B=120, Delta=1 MW at n=3,10,50,100,200. Each uses p>nB and fresh seeded public test masks. The script records inputs, masks, modular wraps, and decoded sums. These masks are not cryptographic randomness.",
             "- W=I-L/(dmax+1). Star and complete spectra are exact; ring/circulant numerical spectra are accompanied by conservative rational bounds derived from sin(pi*t)>=2t for 0<=t<=1/2.",
             "- The fixed public budget uses m*ceil(sqrt(m))*(p-1)*rho_upper^k<1/4. Upward-rounded decimal arithmetic verifies the bound. A conservative IEEE binary64 error bound also keeps the total below 1/2 for every case.",
             "- The arithmetic certificate assumes correctly rounded binary64 addition, multiplication and division, round to nearest, exact small integer inputs/coefficients, and no overflow/underflow. It is not a hardware or compiler formal verification.",
             "- First-correct-decoding rounds are offline oracle diagnostics. They are never used for input-dependent protocol stopping.", "",
             "- Budgets above 12000 updates are evaluated through the known spectral form of W^k, not by executing every update. Other convergence cases execute all updates. Timing runs execute min(k,2000) updates, with their count reported explicitly; no extrapolated full-protocol runtime is claimed. Spectral evaluation is a fast numerical diagnostic, while the separate iteration/arithmetic bounds apply to the actual positive-sum kernel.", "",
             "## Fixed-budget convergence", "",
             "The table gives the maximum final error across three fresh public-mask trials. Directed messages exclude mask exchange, authentication, acknowledgments, and barrier control.", "",
             "| B; Delta MW | Operators | Graph | Max degree | Numerical rho | Certified rho upper | Public rounds | Offline first correct, range | Final scaled error | Directed messages | Timed updates | Median host kernel ms |",
             "|:---|---:|:---|---:|---:|---:|---:|:---|---:|---:|---:|---:|"]
    for c in report["cases"]:
        first = c["convergence"]["offline_first_correct_decoding_rounds"]
        lines.append(f'| {c["input_data"]["B"]}; {c["input_data"]["delta_MW"]} | {c["n_operators"]} | {c["graph"]} | {c["maximum_degree"]} | '
                     f'{c["rho_numerical"]:.7f} | {c["rho_upper_float"]:.7f} | '
                     f'{c["budget"]["public_iterations"]} | {min(first)}–{max(first)} | '
                     f'{max(c["convergence"]["final_max_scaled_errors"]):.5g} | '
                     f'{c["logical_directed_scalar_messages"]:,} | '
                     f'{c["host_kernel_timing"]["measured_updates_per_repeat"]} | '
                     f'{1000*c["host_kernel_timing"]["median_seconds"]:.3f} |')
    lines += ["", "## Reliable loss and delay model", "",
              "Each hypothetical data-plus-ack attempt occupies a 1 or 10 ms slot. Directed attempts fail independently with probability 0, .05, or .20. Successful messages are buffered; only failed messages retry. Every node waits for completion of every message in the round before applying the original synchronous update. A 120-round event-level test at 20% loss verifies every state bitwise against the lossless reference. One private field-element exchange per operator-cycle edge is charged before averaging; channel establishment and cryptographic overhead are excluded.", "",
              "Two schedules expose the hardware assumption: parallel_links permits all links simultaneously; one_peer_per_node uses deterministic greedy edge-color batches, permitting a node to exchange with only one peer at a time. The latter charges the star hub and complete graph for their high degree. It is a specific feasible schedule, not a claim of optimal scheduling.", "",
              "The network Monte Carlo samples the maximum-of-geometric retry law per batch, with a recorded per-episode tail truncation bound below 1e-12. Expected attempts equal delivered scalars/(1-loss). Deadline fractions are finite-sample estimates, with Wilson 95% intervals in JSON; a fraction of 1 is not a probability-one guarantee. The 1, 10, 60, and 300 second deadlines are illustrative coordination budgets; the additional 900-second deadline matches the separate physical study's 15-minute interval. On a missed deadline this protocol makes no output available: it does not decode early. No physical fallback or response to a missed bid is simulated here. The model excludes correlated/burst loss, disconnects, malicious behavior, finite retry limits, clock synchronization, bandwidth/packet-size effects, encryption, and control traffic. Persistent loss can prevent progress; no deadline guarantee follows from the unbounded retry model.", "",
              "Selected fine-grid 20%-loss, 1 ms-slot cases below report expected network time and the Monte Carlo 95th percentile. All 0/.05/.20-loss and 1/10 ms-slot results, including 1-second deadline success, are in the JSON. Host kernel times are kept separate and are not added to synthetic network time.", "",
              "| Operators | Graph | Schedule | Batches/round | Expected seconds | Sample p95 seconds | Fraction <=10 s | Fraction <=60 s | Expected scalar attempts |",
              "|---:|:---|:---|---:|---:|---:|---:|---:|---:|"]
    for c in report["cases"]:
        if c["input_data"]["B"] != 10000:
            continue
        for net in c["network_models"]:
            if net["independent_attempt_loss_probability"] == .2 and net["hypothetical_attempt_slot_ms"] == 1:
                lines.append(f'| {c["n_operators"]} | {c["graph"]} | {net["schedule"]} | '
                             f'{net["batches_per_round"]} | {net["expected_network_seconds"]:.3f} | '
                             f'{net["sample_p95_seconds"]:.3f} | {net["deadline_met_fractions"]["10"]:.3f} | '
                             f'{net["deadline_met_fractions"]["60"]:.3f} | '
                             f'{net["expected_attempted_directed_scalar_transmissions"]:,.1f} |')
    lines += ["", "## Interpretation and manuscript wording", "",
              "Topology changes both the spectral convergence rate and per-round communication. A complete graph reaches the average in one synchronous update but requires quadratic total edges and high per-node degree. A ring has bounded degree but slow convergence. The degree-at-most-four circulant adds long-range links and improves convergence without a growing degree; the utility-centered star retains a degree-n bottleneck. These conclusions concern the specified graph families and gain, not arbitrary networks.", "",
              "Suggested text: 'A reproducible communication study uses 3, 10, 50, 100 and 200 operators and compares star, ring, degree-at-most-four circulant, and complete public consensus graphs. Fixed public iteration budgets combine conservative spectral and binary64 error bounds; all encoded totals are recovered across three fresh public-mask trials per case. Offline decoding thresholds, message counts, and one-host computation times are reported separately. Reliable buffered retransmission preserves the synchronous iteration under independent modeled losses, while changing the completion time. A node-limited matching schedule exposes the star hub degree and distinguishes connectivity from practical communication scalability. Slot duration and deadlines are scenario assumptions, not measured network or microgrid timings.'", "",
              "The public communication graph and private mask graph remain distinct. Improving the former does not improve the latter's collusion tolerance. This study does not extend privacy to constrained bids, metering, active attacks, or dropped participants.", "",
              "Artifacts: `communication_validation.json`, `convergence_n100.svg`, and this report. JSON contains complete inputs, masks, seeds, bounds, sample traces, network assumptions, hardware/software metadata and all results.", ""]
    (directory/"README.md").write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", nargs="+", type=int, default=[3, 10, 50, 100, 200])
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--timing-repeats", type=int, default=3)
    parser.add_argument("--network-episodes", type=int, default=2000)
    parser.add_argument("--skip-fine-grid", action="store_true")
    parser.add_argument("--network-only", action="store_true",
                        help="Refresh synthetic network results while preserving recorded kernel timings")
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).resolve().parents[1]/"output/communication")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.network_only:
        report_path = args.output_dir/"communication_validation.json"
        report = json.loads(report_path.read_text())
        for case in report["cases"]:
            spec = graph(case["n_operators"], case["graph"])
            case["network_models"] = network_suite(spec, case["budget"]["public_iterations"],
                                                   args.network_episodes)
        report_path.write_text(json.dumps(report, indent=2)+"\n")
        write_report(report, args.output_dir)
        write_svg(report["cases"], args.output_dir/"convergence_n100.svg")
        print("Refreshed network results; preserved existing measured kernel timings.")
        return
    cases = []
    profiles = [(n, 120, 1.0) for n in args.sizes]
    if not args.skip_fine_grid:
        profiles += [(n, 10000, .001) for n in (10, 50, 100)]
    for n, bound, delta in profiles:
        assert n >= 3
        for graph_id, kind in enumerate(GRAPHS):
            spec = graph(n, kind)
            data = inputs(n, graph_id, args.trials, bound=bound, delta=delta)
            budget = certified_budget(spec, data["prime"])
            checked = convergence(spec, data, budget)
            timing = kernel_timing(spec, checked.pop("initial"),
                                   budget["public_iterations"], args.timing_repeats)
            nets = network_suite(spec, budget["public_iterations"], args.network_episodes)
            case = {"n_operators": n, "n_consensus_nodes": n+1, "graph": kind,
                    "graph_edges": spec["edges"], "maximum_degree": spec["dmax"],
                    "circulant_shifts": spec["shifts"],
                    "gain": 1/spec["denominator"],
                    "rho_numerical": spec["rho_numerical"],
                    "rho_upper_rational": str(spec["rho_upper"]),
                    "rho_upper_float": float(spec["rho_upper"]),
                    "budget": budget, "input_data": data,
                    "logical_directed_scalar_messages": 2*len(spec["edges"])*budget["public_iterations"],
                    "private_mask_field_elements_per_trial": n,
                    "convergence": checked, "host_kernel_timing": timing,
                    "network_models": nets}
            cases.append(case)
            print(f'B={bound:5d} n={n:3d} {kind:10s}: k={budget["public_iterations"]:6d}, '
                  f'error={max(checked["final_max_scaled_errors"]):.4g}, '
                  f'kernel={timing["median_seconds"]:.4f}s', flush=True)
    retry_check = finite_retry_schedule_check()
    assert any(r["modular_wraps"] > 0 for c in cases for r in c["input_data"]["runs"])
    report = {"schema_version": 1, "status": "all checks passed", "seed": BASE_SEED,
              "environment": {"python": sys.version, "numpy": np.__version__,
                              "platform": platform.platform(), "machine": platform.machine(),
                              "processor": platform.processor()},
              "arguments": {k: str(v) if isinstance(v, Path) else v
                            for k,v in vars(args).items()},
              "scope": {"production_crypto": False, "deployed_network": False,
                        "dynamic_microgrid_physics": False,
                        "arithmetic_certificate": "conditional on standard IEEE binary64 model; see source and README",
                        "mask_security": "public seeded validation data, not secure randomness",
                        "mask_graph": "operator cycle, separate from public consensus graph",
                        "timing": "one-host kernel time and synthetic network time kept separate"},
              "event_level_retry_check": retry_check, "cases": cases}
    (args.output_dir/"communication_validation.json").write_text(json.dumps(report, indent=2)+"\n")
    write_report(report, args.output_dir)
    write_svg(cases, args.output_dir/"convergence_n100.svg")
    print(f'Wrote {args.output_dir}; {len(cases)} topology/size cases, '
          f'{len(cases)*args.trials} fresh-mask trials.', flush=True)


if __name__ == "__main__":
    main()
