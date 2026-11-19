# Draft responses to Reviewer 1, comments 2 and 4

## Comment 2: uncertainty dependence and the price-deviation metric

We have corrected the covariance expansion and added an explicit interpretation of the objective. With net-load forecast error e, aggregate forecast X, and total bid Q, the price difference is 2a2[e − (X − Q)]. Its mean square is therefore 4a2²[Var(e) + (X − Q)²]. The covariance expression is Var(e) = σ_LL − 2Σ_i σ_iL + Σ_ij σ_ij.

This decomposition shows that bids affect the mean price bias, while real-time price variance is independent of day-ahead bids under the assumed affine prices and exogenous real-time generation. Covariance contributes only bid-independent constants to the controller objectives. Hence equilibrium bids are exactly invariant to covariance when forecast means, coefficients, and deterministic bid limits stay fixed. We now state this structural limitation explicitly. The metric measures the common two-price spread entering imbalance settlement, but omits each operator's exposure and loss direction; it is used for tractable objective shaping, and no equivalence to profit variance or tail-loss control is claimed.

We added a four-level uncertainty-scale sweep that reports full expected profits, real-time price variance, and mean-square price deviation. Independent Gaussian quadrature checks the settlement moments and unilateral optimal bids. The expanded version also derives realized-profit variance and reports that ERS increases all three operators' profit standard deviations in the illustrative case. We removed the misleading profit-variance footnote and the implication that the method reduces stochastic price volatility.

## Comment 4: upper bidding limits and constrained Nash/Pareto results

We now distinguish the uncapped model from the game with finite quantity limits 0 ≤ q_i ≤ u_i, including the forecast-relative rule u_i = x_i^o. Theorem 3 proves uniqueness of the constrained Nash equilibrium: q_i* = min(u_i, k_i r), where k_i = 1 + 4a2β_i and r + Σ_i min(u_i, k_i r) = X. The scalar equation accounts for boundary optimality and recomputation of all rival bids, rather than simply clipping an uncapped equilibrium.

The constrained aggregate-profit-maximizing total is Q† = min(X/2, Σ_i u_i). We provide a necessary-and-sufficient condition on the prescribed coefficients for the equilibrium to attain this total and a constructive common-coefficient design. The original uncapped coefficient-sum condition is not imposed on capped equilibria. The old midpoint saturation construction has been replaced; the revised method remains feasible in the five-operator counterexample that gave a negative bid under that construction.

Numerical comparisons distinguish fixed RN, fixed original ERS, and coordinated coefficient retuning. They include cap shortages and zero caps, as well as 200 equilibria across 100 heterogeneous instances checked against unilateral best responses. We clarify that a cap relative to a forecast mean limits commitments but does not guarantee stochastic delivery. The constrained design is an analytical extension and does not inherit the uncapped masked protocol's privacy guarantee.
