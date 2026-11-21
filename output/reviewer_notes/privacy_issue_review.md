# Algorithm 1 privacy issue — review note

Date: November 19, 2026
Status: Documented for later decision. No manuscript changes authorized or made for this note.
Source: main_before_reviewer_changes.tex (original manuscript).

## Reviewer concern

Reviewer 1: “Clarify and rigorously validate the privacy-preserving claim of Algorithm 1; standard consensus may reveal x_i0.”

Related: Reviewer 3, comment 4; Reviewer 4, issue 5.

## What the algorithm computes

The algorithm estimates the operator count and aggregate forecast, then computes each operator's share locally:

- xi_i converges to 1/(n+1).
- psi_i converges to X/(n+1), where X = sum_j x_j^o.
- eta_i = x_i^o xi_i / psi_i converges to x_i^o/X.

Thus the final share is a percentage. This does not mean the exchanged messages contain only percentages.

## Concrete disclosure

Algorithm 1 initializes psi_i(0) = x_i^o and uses neighboring psi states in its consensus updates. Under ordinary direct state exchange, a neighbor receiving the first state receives the individual forecast itself.

Example: x^o = (120,100,80) MW. Operator 1 ultimately computes eta_1 = 40%, but its initial transmitted psi_1 is 120 MW.

The manuscript states that every operator communicates with the utility. If the utility is adjacent to all operators, direct initial-state exchange reveals every forecast to it.

Omitting an initial message does not by itself establish privacy. For a star with utility state psi_{n+1}(0)=0, the first leaf update is psi_i(1)=(1-epsilon_i)x_i^o. If this updated state is observed and epsilon_i is known, the forecast can be reconstructed by division by 1-epsilon_i. This is a reconstruction illustration, not a separate complete protocol for avoiding initial exchange.

## Distinctions to preserve

- Distributed computation: local updates and neighbor-to-neighbor communication. The algorithm can have this property without confidentiality.
- Restricted communication: information is sent only to neighbors, rather than broadcast to every participant. This does not hide it from those neighbors.
- Forecast confidentiality: an observer cannot infer an individual forecast beyond explicitly permitted disclosures. This requires a stated observer/adversary model and a supporting argument.
- Local final calculation: eta_i may remain local, although messages used to obtain its denominator disclose x_i^o.
- Encryption in transit: can protect against external interception, but does not hide message contents from the intended receiving neighbor.

The demonstrated issue concerns confidentiality. It does not itself refute the bidding formulas or establish failure of consensus convergence. Heterogeneous gains and stopping/convergence conditions are separate reviewer issues.

## Claims needing review in the original source

Locate by text rather than relying on line numbers, which may change:

1. Introduction: “implemented in a privacy-preserving and fully distributed manner.”
2. Introduction: “exchanges only privacy-preserving information” and “all market-sensitive data remain local to each operator.”
3. Figure caption: “privacy-preserving consensus variables.”
4. Formulation: “preserving privacy through minimal information exchange.”
5. Opening of distributed section: “privacy is protected.”
6. Algorithm 1 initialization: psi_i(0)=x_i^o.
7. Paragraph after Algorithm 1: “There is no explicit sharing of private information.”
8. Numerical discussion: “without any operator revealing private information.”
9. Broader claims that Nash strategies are preferable for privacy: distinguish information required to determine bids from information exposed by this implementation.

## Options for later decision

### A. Retain ordinary consensus and narrow the claim

Describe the algorithm as a distributed implementation using neighboring state exchange. Explicitly acknowledge that individual forecasts are disclosed to direct recipients and may also be inferable from subsequent observations. Do not claim confidentiality from non-neighbors without analysis either.

This requires targeted wording changes, not a new bidding model. Confidential aggregation against curious participants could be identified as future work. It is appropriate only if forecast confidentiality is not essential to the claimed contribution.

### B. Retain confidentiality as a contribution

Modify the protocol, for example using the previously investigated pairwise masking approach. Specify private channels, permitted outputs, passive/adversarial behavior, collusion/connectivity conditions, repeated executions, encoding, and correctness. Supply the proof supporting the exact guarantee.

This is a substantive protocol addition. Earlier revision experiments exist, but they have NOT been selected for incorporation into the original manuscript under the user's new targeted-edit instructions.

### C. Explicitly trust receiving neighbors

If the intended application genuinely treats all neighbors—including the utility—as trusted recipients, disclose that assumption and describe the actual limited communication property. Trust does not substantiate the original claim that market-sensitive data remain local.

## What additional evidence can and cannot do

A small replay demonstrating initial-message disclosure can illustrate the problem. Code and expanded experiments can be supplied on GitHub later. Neither simulation nor a repository link replaces a privacy definition and proof if a formal privacy guarantee is retained.

## Questions to resolve before editing

- From whom must forecasts be confidential: other operators, the utility, external listeners, or some combination?
- Are neighboring recipients explicitly trusted to know individual forecasts?
- Is privacy essential to the paper, or is limited/local communication the intended contribution?
- Should we correct the wording only, or introduce a protected aggregation protocol?

## Current decision

Document only. Do not modify main_before_reviewer_changes.tex, main.tex, or main_LONG.tex. No final choice between narrowing the claim and changing the protocol has been made.
