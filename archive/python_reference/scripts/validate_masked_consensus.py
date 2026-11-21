#!/usr/bin/env python3
"""Deterministic, executable checks for masked forecast-sum consensus.

This is a validation model, not a cryptographic implementation. random.Random
provides reproducible PUBLIC test masks, not secure production randomness.
The threat model is passive: private authenticated pairwise mask exchange,
fresh independent uniform field masks, and a connected graph induced by the
honest operators after colluders are removed. The utility may observe every
masked initial value and every subsequent consensus message. n is public.

For Delta=1, a public forecast bound B, integer forecasts 0 <= x_i <= B,
and public prime p > n*B, initialize operator i at
y_i = (x_i + sum(outgoing masks) - sum(incoming masks)) mod p and the utility
at zero. Average these residues over a star using gain 1/(n+1). Every node
decodes round((n+1)*psi_i(k)) mod p. The sum S is intentionally disclosed.
No claim is made that output bids themselves reveal nothing about forecasts.

Run from any directory; the default JSON location is relative to this file.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
import math
from pathlib import Path
import random


SEED = 20261117


def masked_values(xs: list[int], edges: list[tuple[int, int]],
                  masks: tuple[int, ...] | list[int], prime: int) -> list[int]:
    values = xs.copy()
    for (source, target), mask in zip(edges, masks, strict=True):
        values[source] += mask
        values[target] -= mask
    return [value % prime for value in values]


def is_prime(value: int) -> bool:
    return value >= 2 and all(value % d for d in range(2, math.isqrt(value) + 1))


def next_prime(lower: int) -> int:
    while not is_prime(lower):
        lower += 1
    return lower


def round_fraction(value: Fraction) -> int:
    """Nearest integer for nonnegative exact values; no tie occurs in checks."""
    assert value >= 0
    return (2 * value.numerator + value.denominator) // (2 * value.denominator)


def iteration_budget(n: int, prime: int) -> tuple[int, float]:
    """Certify m^(3/2)(p-1)rho^k < 1/4 by an integer inequality.

    For z(0) in [0,p-1]^m, ||z(0)-mean(z(0))*1||_2 <= sqrt(m)(p-1).
    The error of m*z_i(k) is consequently at most
    m^(3/2)(p-1)*rho^k, rho=n/(n+1). Squaring makes the certification exact.
    This is an ideal-arithmetic convergence bound, not a floating-point bound.
    """
    m = n + 1
    rho = n / m
    constant = m ** 1.5 * (prime - 1)
    k = max(1, math.ceil(math.log(0.25 / constant) / math.log(rho)))
    while 16 * m**3 * (prime - 1)**2 * n**(2*k) >= m**(2*k):
        k += 1
    return k, constant * rho**k


def star_consensus(values: list[int], k: int, exact: bool) -> list:
    """Synchronous star updates: leaves first and utility last."""
    n = len(values)
    m = n + 1
    if exact:
        z = [Fraction(x) for x in values] + [Fraction(0)]
        gain = Fraction(1, m)
    else:
        z = [float(x) for x in values] + [0.0]
        gain = 1.0 / m
    for _ in range(k):
        leaves, hub = z[:-1], z[-1]
        z = [(1-gain)*leaf + gain*hub for leaf in leaves] + [
            (1-n*gain)*hub + gain*sum(leaves)
        ]
    return z


def equilibrium_from_sum(xs: list[int], total: int) -> dict:
    """Interior equal-risk-seeking profile of the paper, at a2=0.1."""
    n = len(xs)
    beta = Fraction(1-n, 4*n) / Fraction(1, 10)
    q = Fraction(total, 2*n)
    alpha = [q / x for x in xs]
    return {
        "beta": float(beta), "bid_mw_per_operator": float(q),
        "alpha": [float(value) for value in alpha],
        "all_bids_within_forecast": all(0 <= value <= 1 for value in alpha),
    }


def correctness_case(n: int, seed: int, exact: bool = False) -> dict:
    rng = random.Random(seed)
    public_bound = 120
    xs = [120, 100, 80] if n == 3 else [rng.randint(80, 120) for _ in range(n)]
    total = sum(xs)
    prime = 1009 if n == 3 else next_prime(n * public_bound + 1)
    assert all(0 <= value <= public_bound for value in xs)
    assert prime > n * public_bound
    # A cycle of private mask links; for n=3 it is a triangle. The PUBLIC
    # averaging graph is independently a star through the observing utility.
    edges = [(i, (i+1) % n) for i in range(n)]
    masks = [rng.randrange(prime) for _ in edges]
    values = masked_values(xs, edges, masks, prime)
    assert sum(values) % prime == total
    k, bound = iteration_budget(n, prime)
    z = star_consensus(values, k, exact=exact)
    residue_sum = sum(values)
    scaled = [(n+1)*value for value in z]
    rounded = [round_fraction(value) if exact else round(value) for value in scaled]
    decoded = [value % prime for value in rounded]
    assert all(value == residue_sum for value in rounded)
    assert all(value == total for value in decoded)
    central = equilibrium_from_sum(xs, total)
    distributed = equilibrium_from_sum(xs, decoded[0])
    assert central == distributed
    assert central["all_bids_within_forecast"]
    max_error = max(abs(value-residue_sum) for value in scaled)
    if exact:
        assert max_error < Fraction(1, 4)
        exact_error = {"numerator": str(max_error.numerator),
                       "denominator": str(max_error.denominator)}
    else:
        exact_error = None
    return {
        "n_operators": n, "n_consensus_nodes": n+1, "seed": seed,
        "delta_mw": 1, "public_forecast_bound_B": public_bound,
        "prime": prime, "prime_exceeds_nB": True, "input_forecasts_mw": xs,
        "sum_forecasts_mw": total, "mask_graph_edges": edges,
        "public_test_masks": masks, "masked_initial_values": values,
        "sum_residues_before_modulo": residue_sum,
        "mask_sum_cancels_modulo_prime": True,
        "gain": 1/(n+1), "rho": n/(n+1), "iterations": k,
        "ideal_arithmetic_total_error_bound": bound,
        "ideal_bound_below_quarter_certified_by_integer_inequality": True,
        "arithmetic": "exact rational" if exact else "binary64 empirical check",
        "max_scaled_consensus_error": float(max_error),
        "exact_error_fraction": exact_error,
        "decoded_sums_at_all_nodes": decoded,
        "all_nodes_decode_exact_sum": True,
        "centralized_ers": central, "distributed_ers": distributed,
        "ers_profile_matches_exactly": True,
    }


def view_distribution(xs: list[int], edges: list[tuple[int, int]], prime: int,
                      colluders: tuple[int, ...] = ()) -> Counter:
    """Full initial public messages plus colluders' inputs and incident masks.

    All future deterministic consensus transcripts are functions of this view,
    so identical distributions here imply identical full-consensus transcripts.
    """
    distribution = Counter()
    for masks in product(range(prime), repeat=len(edges)):
        known_masks = tuple(mask for edge, mask in zip(edges, masks, strict=True)
                            if any(node in colluders for node in edge))
        view = (tuple(xs[i] for i in colluders), known_masks,
                tuple(masked_values(xs, edges, masks, prime)))
        distribution[view] += 1
    return distribution


def privacy_exhaustive_checks() -> dict:
    prime = 7
    public_bound = 2
    a, b = [1, 1, 2], [1, 2, 1]
    assert prime > len(a) * public_bound
    assert all(0 <= value <= public_bound for value in a + b)
    triangle = [(0, 1), (0, 2), (1, 2)]
    utility_a = view_distribution(a, triangle, prime)
    utility_b = view_distribution(b, triangle, prime)
    assert utility_a == utility_b
    assert len(utility_a) == prime**2
    assert set(utility_a.values()) == {prime}
    collusion_a = view_distribution(a, triangle, prime, (0,))
    collusion_b = view_distribution(b, triangle, prime, (0,))
    assert collusion_a == collusion_b

    disconnected = [(0, 1)]
    disconnected_a = view_distribution(a, disconnected, prime)
    disconnected_b = view_distribution(b, disconnected, prime)
    assert not (disconnected_a.keys() & disconnected_b.keys())
    components_a = [(a[0]+a[1]) % prime, a[2] % prime]
    components_b = [(b[0]+b[1]) % prime, b[2] % prime]
    for distribution, expected in [(disconnected_a, components_a),
                                   (disconnected_b, components_b)]:
        for _, _, values in distribution:
            assert [(values[0]+values[1]) % prime, values[2]] == expected

    hub_mask_graph = [(0, 1), (0, 2)]
    broken_a = view_distribution(a, hub_mask_graph, prime, (0,))
    broken_b = view_distribution(b, hub_mask_graph, prime, (0,))
    assert not (broken_a.keys() & broken_b.keys())
    # With operator 0 colluding, it knows both masks and reconstructs both
    # singleton honest components: x1=(y1+r01) mod p, x2=(y2+r02) mod p.
    for distribution, xs in [(broken_a, a), (broken_b, b)]:
        for _, known_masks, values in distribution:
            assert [(values[1]+known_masks[0]) % prime,
                    (values[2]+known_masks[1]) % prime] == xs[1:]
    return {
        "prime": prime, "public_forecast_bound_B": public_bound,
        "prime_exceeds_nB": True, "input_a": a, "input_b": b, "common_sum": sum(a),
        "triangle_utility_observer": {
            "enumerated_mask_assignments_per_input": prime**3,
            "distinct_views": len(utility_a),
            "multiplicity_per_view": prime,
            "view_distributions_identical": True,
        },
        "triangle_utility_and_operator_0_collude": {
            "honest_induced_graph_edges": [[1, 2]],
            "enumerated_mask_assignments_per_input": prime**3,
            "distinct_full_views": len(collusion_a),
            "view_distributions_identical": True,
        },
        "disconnected_honest_graph": {
            "edges": disconnected, "components": [[0, 1], [2]],
            "component_sums_input_a": components_a,
            "component_sums_input_b": components_b,
            "view_supports_disjoint": True,
            "component_sums_recovered_for_every_mask_assignment": True,
        },
        "operator_0_collusion_disconnects_honest_graph": {
            "edges": hub_mask_graph, "colluders": [0],
            "honest_components": [[1], [2]],
            "view_supports_disjoint": True,
            "both_honest_inputs_recovered_for_every_mask_assignment": True,
        },
    }


def plain_consensus_attack() -> dict:
    xs = [120, 100, 80]
    gain = Fraction(1, 4)
    first_step = [(1-gain)*value for value in xs]
    recovered = [round_fraction(value/(1-gain)) for value in first_step]
    assert recovered == xs
    plain = star_consensus(xs, 37, exact=True)
    plain_decoded = [round_fraction(4*value) for value in plain]
    assert plain_decoded == [sum(xs)] * 4
    return {
        "initial_messages": xs,
        "initial_messages_reveal_forecasts_exactly": True,
        "utility_initial_value": 0,
        "first_leaf_messages_if_initial_round_hidden": [float(x) for x in first_step],
        "forecasts_reconstructed_from_first_leaf_messages": recovered,
        "plain_consensus_decoded_sum": plain_decoded,
        "centralized_ers": equilibrium_from_sum(xs, sum(xs)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1]
                        / "output/privacy/masked_consensus_validation.json")
    args = parser.parse_args()
    cases = [correctness_case(n, SEED+n, exact=(n == 3)) for n in (3, 10, 50, 100)]
    wraparound = correctness_case(3, 20261121, exact=True)
    assert wraparound["sum_residues_before_modulo"] == (
        wraparound["sum_forecasts_mw"] + 2 * wraparound["prime"])
    report = {
        "schema_version": 1,
        "status": "all deterministic checks passed",
        "reproducible_seed_base": SEED,
        "scope": {
            "protocol": "finite-field pairwise masks; real star averaging of residues",
            "known_n": True, "delta_mw": 1,
            "exact_recovery_condition": "public B bounds every integer forecast; public prime p > n*B guarantees 0 <= S < p",
            "prime_selection": "from public n and B, never from the private realized sum",
            "mask_security": "fresh independent uniform masks shared privately between endpoints",
            "privacy": "passive utility/operator collusion; connected honest mask graph; sum disclosed",
            "not_tested": ["active message tampering", "secure channel implementation",
                           "production cryptographic randomness", "delays or packet losses",
                           "output bid leakage", "performance benchmark"],
            "floating_point_limit": "n=3 uses exact Fraction arithmetic; larger n are empirical binary64 checks, not certified finite-precision implementations",
            "quantization_limit": "Delta=1 inputs here are integers; general rounding introduces up to n*Delta/2 aggregate quantization error",
            "evidence_limit": "finite exhaustive examples and deterministic correctness checks are not a general privacy proof",
        },
        "unmasked_forecast_attack": plain_consensus_attack(),
        "correctness_cases": cases,
        "exact_rational_wraparound_case": wraparound,
        "small_field_privacy_checks": privacy_exhaustive_checks(),
    }
    assert cases[0]["centralized_ers"] == report["unmasked_forecast_attack"]["centralized_ers"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"All checks passed; wrote {args.output}")
    for case in cases:
        print(f"n={case['n_operators']}: p={case['prime']}, k={case['iterations']}, "
              f"S={case['sum_forecasts_mw']}, max scaled error={case['max_scaled_consensus_error']:.6g}")
    print(f"Exact wraparound: residues sum={wraparound['sum_residues_before_modulo']} "
          f"= S+2p, decoded S={wraparound['sum_forecasts_mw']}")


if __name__ == "__main__":
    main()
