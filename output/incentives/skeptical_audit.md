# Independent audit of the proposed settlement extension

Scope: current local main.tex and main_LONG.tex, fixed forecasts, exogenous real-time output, affine prices, and unilateral quantity choices. This is a mathematical audit, not a legal assessment. No manuscripts changed.

## 1. Fixed-input implementation is sound

For n>=2, fix a feasible target t, T=sum(t), and c_i=2a2(X-T-t_i). Let

    transfer_i=-c_i(q_i-t_i)+sum_{j!=i} c_j(q_j-t_j)/(n-1).

Transfers sum to zero for every bid vector, not merely at equilibrium, and vanish individually at q=t. Actual expected profit plus transfer has derivative 2a2(X-Q-q_i)-c_i and own second derivative -4a2. This is an exact potential game with potential

    Phi(q)=sum_i (2a2 X-c_i)q_i-a2 Q^2-a2 sum_i q_i^2.

Its Hessian -2a2(I+11') is negative definite. The target has zero gradient and is the unique Nash equilibrium on any product of closed finite bid intervals containing it. This includes boundary targets. More generally, a one-sided KKT construction could avoid unnecessary charges on binding caps, but is a different transfer choice.

Thus an enforceable settlement turns the assigned objective into actual monetary bidding incentives. It does not require owners to possess negative risk preferences. The selected controller coefficients can remain a design route to the target, but transfers rather than intrinsic risk preferences implement owner incentives. The mechanism can implement arbitrary feasible targets, so its novelty must not be overstated.

## 2. Enforceability is not eliminated

The result assumes the settlement rule cannot be refused after bidding, actual bids are observable, payments are collectible, and identity/membership cannot be changed to avoid charges. A binding aggregator agreement with delegated settlement or collateral is a concrete setting. Cryptographic consensus alone supplies none of these properties. This is conditional incentive-compatible bidding, not enforcement-free implementation.

For bounded bid boxes, worst-case charges can be bounded in advance and collateralized; no such bound exists for an uncapped strategy space. This practical fact favors stating the implementation theorem for capped bids. Liquidity, bankruptcy, metering disputes, and legal authority are not modeled.

Participation must be tested at a separate enrollment stage against a specified opt-out continuation. The zero-on-target transfers preserve the original ERS allocation and profits, but cannot improve its participation comparison merely by being present.

## 3. Coalition robustness: symmetric versus capped

For symmetric uncapped ERS, t_i=X/(2n), c_i=a2 X(n-1)/n, and delta_i=q_i-t_i, direct expansion gives

    U_i(q)-U_i(t)=-2a2 delta_i sum_j delta_j.

Hence if only coalition S changes bids, its total gain is -2a2(sum_{i in S}delta_i)^2 <=0. No coalition can make all members weakly better off and one strictly better off through these fixed-input quantity changes, even with internal budget-balanced side payments. Sum-preserving changes leave all coalition members indifferent. This is a limited coalition result, not a claim covering reporting, exclusion, renegotiation, or changing the external market.

An initial informal coalition calculation was incorrect because it omitted cross-bid effects on market revenue; the corrected direct calculations are in skeptical_checks.py/json.

The heterogeneous capped extension need not share this property. In the manuscript's capped case, X=250, a2=.1, target=(47.5,47.5,30), caps=(120,100,30), and charges=(15.5,15.5,19). Coalition {1,3} changes only q1 to 51.875. Operator 1's net payoff falls by 3.828125 dollars/h, operator 3's rises by 7.65625 dollars/h, and their aggregate rises by 3.828125 dollars/h. An internal side payment can make both better off. All bids remain feasible. Therefore the capped theorem should claim unique unilateral Nash implementation, not general coalition-proofness.

## 4. Information and privacy requirements

The theorem fixes X, caps, and targets before bidding. Strategic reporting is an additional stage; encrypted or masked aggregation protects confidentiality, not truthfulness, and the original passive-adversary privacy guarantee does not validate malicious inputs.

Do not automatically claim profitable forecast lying in the symmetric uncapped case. If the only manipulable input is aggregate forecast H, the mechanism uses t_i=H/(2n), all common charges c=a2 H(n-1)/n, and the subsequent bidding equilibrium is interior, all bids equal [X-H(n-1)/(2n)]/(n+1). Transfers then cancel individually. Every owner's bid-dependent payoff along these symmetric equilibria is maximized when H=X. Thus truthful aggregate reporting is weakly optimal in this restricted report-and-re-equilibrate game; individual reports are not identified because any report vector with correct total is equivalent. Caps, altered membership, authentication, constraints on reports, and other settlement obligations can change this conclusion. It is not a general strategyproofness proof.

For capped targets, private caps are required for target design; the existing masked protocol only handles uncapped aggregate forecasts. Do not extend its privacy theorem to cap-dependent target computation without a new protocol. Actual bids are already submitted to the market, but whether an aggregator can observe every bid for transfers must be specified. A distributed payment computation is also not guaranteed merely by consensus convergence.

## 5. Consumer and efficiency objection remains material

Budget-balanced internal transfers do not change external prices or consumer payments at the target. They implement producer-profit coordination rather than establish social efficiency.

If consumers settle forecast demand L0 day-ahead and deviations in real time, their expected payment is lambda_DA L0+E[lambda_RT(L-L0)]. Under this manuscript's exogenous real-time model, reducing aggregate bids from RN nX/(n+1) to ERS X/2 raises their expected payment by

    a2 L0 X(n-1)/(n+1).

In baseline n=3,X=300,L0=400,a2=.1 this is 6000 dollars/h. Renewable aggregate profit rises 1125 dollars/h. Under consistent two-settlement conventional revenues with physical generation cost independent of bids, the remaining 4875 dollars/h goes to conventional profit. This is an accounting example conditional on that buyer/conventional settlement specification, not a result already stated in the manuscript. Inelastic demand and bid-independent real-time physical dispatch leave physical surplus unchanged; higher producer profits are redistribution, not demonstrated welfare improvement. Consumer utility, elastic demand, and realistic dispatch may produce different welfare conclusions.

This makes Reviewer 3's withholding concern stronger rather than removing it. The mechanism should be framed as a contracted-aggregator implementation of a producer objective, accompanied by transparent consumer/welfare analysis. No legal characterization is made here.

## Recommendation

Add a bounded-domain settlement theorem only if the paper explicitly defines binding settlement authority, observable bids, fixed/verified design inputs, and participation timing. Keep voluntary opt-out and any broader coalition/privacy claims separate. Preserve the existing ERS target theorem as the design result; describe settlement as an implementation layer. A consumer-payment accounting section should accompany it rather than treating budget balance as social desirability.
