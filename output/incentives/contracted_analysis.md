# Contracted-aggregator incentive analysis

Scope: current affine-price expected-profit model, Pi_i=C_i+2a q_i(X−Q), a=a_2>0, fixed authenticated X>0, n≥2, independent box constraints 0≤q_i≤u_i. C_i is bid-independent. Physical feasibility must cover the whole box as already required by the manuscripts. This analysis does not alter either manuscript.

## 1. Implementation of the uncapped ERS target

Let q*=X/(2n), c=aX(n−1)/n, and
T_i(q)=−c(q_i−q*)+[c/(n−1)] Σ_{j≠i}(q_j−q*).
Owners maximize U_i=Pi_i+T_i; beta is no longer their actual utility parameter. The settlement implements the controller's target, rather than making the owner's intrinsic preferences equal to the controller objective.

Budget balance: summing T_i gives zero for every q. Transfers vanish at the target. If δ_i=q_i−q* and D=Σδ_i, direct expansion gives the particularly useful identity
U_i(q)−U_i(q*1)=−2a δ_i D.
Thus any unilateral quantity deviation against target bids loses exactly 2aδ_i². In Scenario 1, a=.1, X=300, n=3, c=20 dollars/MWh: moving 50→100 MW gains $500/h market revenue, pays $1000/h, and loses $500/h net.

Uniqueness: the game has exact potential
Phi(q)=Σ_i(2aX−c)q_i−aQ²−aΣ_i q_i².
Its Hessian is −2a(I+11^T), negative definite. Every Nash equilibrium is the unique maximizer of this strictly concave potential over the product feasible set. If the ERS target is feasible it is therefore the unique Nash equilibrium, including on a finite physical bid box. This is a Nash guarantee, not a dominant-strategy guarantee: best responses generally depend on others' bids.

A strategic beta setting that merely changes an operator's submitted quantity cannot help against the target bids. Full equivalence to arbitrary dynamic controller games is not automatic. For the uncapped symmetric case, however, committing to q_i=q*+δ and allowing all other owners to re-equilibrate under this same settlement gives each rival q*=X/(2n)−δ/n (while interior). The deviator's net gain is −2aδ²/n, also nonpositive. This calculation depends on rivals maximizing net monetary payoff, not their original J objectives.

## 2. Coalitions: a useful limited property

For any coalition S changing quantities while outsiders remain at the target, Σ_{i∈S}[U_i(q)−U_i(q*)]=−2a(Σ_{i∈S}δ_i)²≤0. Hence even transferable-utility coalition deviations cannot increase the coalition's total net payoff in this fixed-membership settlement game. If coalition total is unchanged, every member's payoff is unchanged. This is not a claim about coalitions withdrawing, manipulating X, altering the settlement agreement, or changing market institutions.

## 3. Capped target generalization

Take any feasible efficient target q† with Q†=min(X/2,Σu_i). Define c_i=2a(X−Q†−q_i†) and
T_i=−c_i(q_i−q_i†)+(1/(n−1))Σ_{j≠i}c_j(q_j−q_j†).
Again transfers sum to zero for every profile and vanish at target. Net own derivative is 2a(X−Q−q_i)−c_i, zero at q†. The same strictly concave potential with individual c_i proves unique Nash implementation on the bid box, even if target coordinates are at bounds. Because q_i†≤Q†≤X/2, c_i≥0. Against target rivals the exact unilateral gain is −2aδ_i². No claim of coalition proof is made for heterogeneous c_i.

A simpler common-charge version implements the manuscript's water-filled target q_i†=min(u_i,t): choose c=2a(X−Q†−t). At free coordinates the derivative is zero; at upper-bound coordinates it is 2a(t−u_i)≥0. The same potential proves uniqueness. For U≤X/2, choose t=max u_i; RN already attains the efficient all-cap profile, so adding a mechanism is unnecessary. These common charges require X, t and the aggregate target; arbitrary individual-target charges require additional target information. The existing private consensus protocol does not privately compute the constrained target.

## 4. Voluntary entry with unanimous commitment and RN fallback

An explicit two-stage arrangement can support voluntary precommitment:
1. Authenticated means/caps and the settlement are fixed before bidding.
2. All n must consent. If anyone rejects, every operator plays the ordinary RN game and no transfers apply.
3. After unanimous consent, membership and settlement are binding for that interval; owners remain free to choose quantities but cannot avoid settlement.

In the uncapped case, the existing ERS−RN gain is aX²(n−1)²/[2n(n+1)²]>0 for every owner. Since transfers vanish at ERS, unanimous consent is a strict best response when all others consent. This gives a subgame-perfect equilibrium with unanimous entry and subsequent unique ERS bidding under credible RN continuation. It does not yield unique entry: if another operator rejects, one's own entry decision may have no effect. For n≥2, coordination-failure rejection equilibria can exist. Nor does this establish ex post unrestricted opt-out stability. The precise outside option is RN for everyone, not one outsider facing an operating ERS coalition.

Capped ERS need not improve every owner's market payoff. Let p_i† and p_i^RN be expected market profits under the chosen target and RN equilibrium; S=Σ_i(p_i†−p_i^RN)≥0. A bid-independent fixed transfer
b_i=p_i^RN−p_i†+w_i S, with w_i≥0 and Σw_i=1,
is budget balanced and leaves bidding incentives unchanged. At target, net payoff becomes p_i^RN+w_i S. With S>0 and all w_i>0, every owner strictly prefers unanimous entry to the specified RN fallback. If S=0, only weak participation is obtained. Bid-independent C_i cancels from these differences, so forecast covariances do not need to be revealed merely to calculate b_i in this model. Computing capped RN and target quantities still needs cap information or an additional secure computation protocol.

## 5. Essential limits for the response to reviewers

- An enforceable settlement or escrow is still required. No mechanism makes monetary obligations binding by mathematical fiat. This replaces obedience to a controller with self-interested optimal bidding inside a binding contract; it is not 'without enforcement'.
- Uniform treatment of acceptance and fallback is a specified institutional arrangement, not something derived from the electricity market model. Without it, participation must be recomputed.
- X, caps, membership, and targets are fixed/authenticated. Privacy does not imply truthful reporting. A forecast-reporting stage is an additional game and not solved here.
- Strict individual entry guarantees presume zero/negligible participation costs; costs must be deducted from surplus before claiming individual rationality.
- Budget balance is among operators, excludes administrative/enforcement costs, and provides no consumer-welfare guarantee. The market-power/withholding objection remains separate.
- Do not pay both the original artificial risk objective and this transfer as if both were actual owner preferences: the claim is implementation of its target through actual monetary incentives.
- Forecast means and physical trajectories must remain exogenous to bidding as in the manuscript. The formula does not extend automatically to bid-following storage or binding conventional dispatch.

Recommendation: add a scoped proposition on budget-balanced implementation and a separately stated unanimous-entry corollary with RN fallback. Describe it as 'incentive-compatible bidding under an enforceable aggregator settlement', not unrestricted self-enforcement. Leave consumer welfare, truthful reports, and general opt-out behavior explicit.

Independent-audit capped counterexample (cross-checked algebraically): X=250, a=.1, q†=(47.5,47.5,30), c=(15.5,15.5,19), equal off-diagonal rebate shares. Coalition {1,3} increases only q1 by 4.375 MW. Agent 1 loses 3.828125 dollars/h, agent 3 gains 7.65625 dollars/h, giving coalition total gain 3.828125 dollars/h. The cap q1≤120 does not obstruct this move. Thus the heterogeneous-charge capped mechanism is not coalition-proof with side payments, despite unique unilateral Nash implementation. The symmetric uncapped coalition identity must not be extended to this case.
