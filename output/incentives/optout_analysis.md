# Opt-out audit

Model checked against main.tex: Pi_i=C_i+2a qi(X-Q), a=a2>0. X and caps are fixed verified inputs; C_i does not change with bids or transfers. Target q*=X/(2n), target net profit C_i+aX²/(2n), and proposed transfers vanish there. No claim below covers strategic reports, coalition deviations, or alternative market clearing.

## 1. Exit while rivals retain target quantities

An exiting owner owes no transfer. Its unconstrained best response is b=X(n+1)/(4n), exceeding its target for n>1. Its gain is aX²(n-1)²/(8n²)>0. With cap u_i>=q*, let delta=min(u_i,b)-q*. The gain is aX(n-1)delta/n-2a delta², strictly positive whenever u_i>q*. Thus feasibility of ERS does not imply participation.

Baseline n=3, X=300, a=.1, caps=(120,100,80): target=(50,50,50), variable profit=1500 per owner. Operator 1 exits, bids100, earns2000: gain500 dollars/hour. Operator3 can bid80 and gains420. Constants cancel.

## 2. Exit and contractual insiders re-equilibrate

Specify continuation carefully: retain insiders' linear own-bid charge c=aX(n-1)/n, and redistribute their charges only among remaining insiders. With at least two insiders this is budget-balanced and leaves their own marginal incentive unchanged. The exiting owner neither pays nor receives these transfers. This is an explicit continuation rule, not an automatic consequence of the original n-player formula.

Write d=c/(2a)=X(n-1)/(2n). With no binding caps, outsiders best respond q0=r=X-Q and insiders qj=r-d. Hence
r=X(n²+1)/(2n(n+1)), qj=X/[n(n+1)].
The outsider's variable profit is 2ar²; compared to aX²/(2n), the gain is positive for n>1 because (n²+1)²-n(n+1)²=n⁴-n³-n+1=(n-1)²(n²+n+1)>0.

Baseline unconstrained continuation: (125,25,25), outsider variable profit3125, gain1625. With operator1 cap120 and other caps100,80, equilibrium=(120,80/3,80/3). Outsider profit3040, gain1540. These are different from the paper's beta-deviation continuation, since remaining owners here have fixed linear settlement charges, not unchanged beta objectives.

For n=2 a lone remaining insider cannot fund its own nonzero incentive through internal budget-balanced redistribution. A rule retaining that charge needs an external counterparty or revised mechanism. Do not apply the insiders-only construction blindly to n=2.

## 3. Why static participation transfers cannot fix this outside option for everyone

Under the fixed-rivals continuation, all uncapped owners have strictly positive exit gains. Keeping the original quantity target and requiring each owner prefer participation would require target transfers T_i>=gain_i>0. Their sum is positive, contradicting budget balance. In particular, zero target transfers cannot work. This is a conditional impossibility for this outside option and target, not an impossibility theorem for all contracts or repeated games. Binding caps can eliminate some gains; if every target equals its binding cap, this argument does not establish failure.

For a general constrained target qdag, fixed-rivals outside best response is clip((X-Qdag+qdag_i)/2,0,u_i). Its exact profit advantage is 2a[b(X-Qdag+qdag_i-b)-qdag_i(X-Qdag)]. This expression provides a participation test independent of the target implementation mechanism.

## 4. Unanimous, all-or-nothing precommitment

A different feasible institutional rule: all operators must commit before bidding; one rejection cancels the whole arrangement and everyone plays the RN bidding game. With enforceable interval commitments, fixed verified inputs, uncapped ERS feasible and zero transfers at target, every owner prefers unanimous ERS to all-RN by
 aX²(n-1)²/[2n(n+1)²]>0.
Baseline gain375 dollars/hour per owner. Thus unanimous acceptance is a strict unilateral entry equilibrium, followed by incentive-compatible contractual bidding. Failure-to-form equilibria can also exist: if others reject, one owner's acceptance alone may change nothing. There is no unique entry equilibrium established here.

This mechanism changes the outside option: refusal destroys the arrangement instead of allowing free riding on participating rivals. Credibility requires a common launch gate and enforceable precommitment (or delegated bidding and settlement authority). It is not enforcement-free self-enforcement. Entry costs subtract from gains; capped heterogeneous targets require separate individual participation checks, since aggregate efficiency alone gives no individual guarantee.

## Recommendation

Use the contracted aggregator plus an explicit unanimous launch/fallback rule if institutionally appropriate. State conditional bidding Nash implementation separately from the participation equilibrium. Reject the claim that the current linear settlement alone supports unrestricted opt-out. Do not claim voluntary participation merely from ERS versus RN comparisons unless all-RN is the specified credible fallback.
