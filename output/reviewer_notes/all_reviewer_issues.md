# Consolidated reviewer issue register

Updated: November 19, 2026. Purpose: issues and response options for later review, not a final rebuttal.

## Working instructions

- Working original: `main_before_reviewer_changes.tex`.
- Do not modify `main_LONG.tex`.
- Preserve the original manuscript and contribution; make only approved, targeted changes.
- Decide whether each concern needs a correction, scope clarification, future work, or supplementary evidence.
- Earlier revised manuscripts and analyses are reference material, not approved wholesale replacements.
- No manuscript edits are made by creating this register. The privacy and uncertainty issues have been discussed against the original; no response option has yet been selected for those issues.

## All reviewer comments and issues

### Reviewer 1

#### R1.1

**Review comment:** Clarify and rigorously validate the privacy-preserving claim of Algorithm 1; standard consensus may reveal x_i0

**Issue:** Ordinary consensus is distributed but does not conceal forecasts from recipients. Algorithm 1 initializes psi_i(0)=x_i^o; a neighbor receiving that state receives the forecast. Computing a percentage locally at the end does not undo disclosure. On a star, a known-gain first updated leaf state can also reveal the forecast.

**Response options / remaining work:** Required claim correction or protocol change. Decide whom data must be hidden from. Either describe limited neighboring communication without a confidentiality guarantee, or add a protected aggregation protocol and proof. A GitHub replay can illustrate disclosure but cannot replace the argument.

#### R1.2

**Review comment:** Explain why the proposed “risk-aware” strategy is largely independent of uncertainty magnitude, and justify the price-deviation risk metric.

**Issue:** The price-deviation metric decomposes into real-time price variance plus squared mean-price bias. At fixed means, coefficients, and bid limits, uncertainty contributes bid-independent constants to both expected profit and the objective. Hence optimal bids do not respond to covariance magnitude; the same rules arise in the deterministic case. The metric ignores individual imbalance exposure and loss direction and can be positive without uncertainty.

**Response options / remaining work:** Required clarification of interpretation and claim scope. Justify a quadratic settlement-price-deviation objective, not comprehensive downside risk protection. Changing to an uncertainty-responsive risk objective would change the model and may belong in future work. Sensitivity experiments can supplement the explanation on GitHub.

#### R1.3

**Review comment:** Clarify whether beta_i is an exogenous risk preference or an endogenous strategic/design variable, and analyze unilateral deviations in beta_i.

**Issue:** The original text mixes intrinsic owner preferences with coefficients chosen to coordinate controllers. Selecting a controller objective does not establish that a profit-maximizing owner will retain it. Coefficient changes with re-equilibrating rivals differ from direct quantity deviations.

**Response options / remaining work:** Required interpretation decision and targeted analysis. Prior revisions explored prescribed controller coefficients and profitable deviations. They are reference material, not automatically approved edits to the original.

#### R1.4

**Review comment:** Resolve the inconsistency regarding alpha_i>1 and formally incorporate practical bidding limits into the Nash/Pareto analysis.

**Issue:** Uncapped equilibrium formulas can produce alpha_i>1. A finite bid constraint changes best responses and potentially the equilibrium and efficient coefficient design; clipping an uncapped solution is insufficient.

**Response options / remaining work:** Required feasible-set clarification and analysis if practical capped applicability is claimed. Earlier revisions contain constrained equilibrium and retuning results. Review for selective incorporation. Large random checks may be supplementary.

#### R1.5

**Review comment:** Revisit the scalability claim because the consensus gain and star-graph node degree depend on n.

**Issue:** Consensus gain restrictions depend on maximum degree. In a star, fixed algebraic connectivity alone does not imply size-independent convergence: with common h=1/(n+1), the contraction factor is n/(n+1). Hub communication and scheduling also grow with size.

**Response options / remaining work:** Required correction to scalability wording and bounds. Earlier graph/message/deadline studies can support a concise manuscript result, with detailed evidence on GitHub.

#### R1.6

**Review comment:** Strengthen the connection to dynamic microgrids and cyber-physical microgrid scope; the current formulation is primarily an electricity-market model.

**Issue:** A static electricity-market game is not a model of electrical microgrid dynamics. Battery energy evolution demonstrates supervisory operation but does not establish voltage/frequency stability or network feasibility.

**Response options / remaining work:** Clarify application scope. Earlier synthetic sequential battery studies offer supporting evidence, but network power flow, fast dynamics, and strategic storage remain separate extensions unless explicitly added.

#### R1.7

**Review comment:** Provide substantially stronger numerical validation, including larger systems, uncertainty sensitivity, and communication/convergence studies.

**Issue:** A three-operator example principally illustrates closed-form results; it does not demonstrate large-system implementation or practical network behavior.

**Response options / remaining work:** Supporting evidence needed. Earlier work includes larger synthetic systems, uncertainty/cap checks, and communication studies. Preserve the distinction between constrained physical simulations and uncapped private consensus. Put essential results in the manuscript and full reproducible runs on GitHub; no repository publication has been performed for this note.

### Reviewer 2

#### R2.1

**Review comment:** The role of the risk parameter is not clear. The paper treats $\beta_i$ both as a generator's risk preference and as a design variable that can be adjusted to obtain a desirable equilibrium. These are different interpretations. If it represents a genuine preference, the market operator cannot simply choose it.

**Issue:** The original text mixes intrinsic owner preferences with coefficients chosen to coordinate controllers. Selecting a controller objective does not establish that a profit-maximizing owner will retain it. Coefficient changes with re-equilibrating rivals differ from direct quantity deviations.

**Response options / remaining work:** Required interpretation decision and targeted analysis. Prior revisions explored prescribed controller coefficients and profitable deviations. They are reference material, not automatically approved edits to the original.

#### R2.2

**Review comment:** Showing that every generator earns more when all adopt ERS does not show that unilateral deviation is unprofitable. Using Table I, ERS gives day-ahead quantities (50,50,50) MW. A generator maximizing expected profit can increase its own quantity from 50 to 100 MW, with the other two unchanged, and gain another $500/h.

**Issue:** The reviewer is correct: at X=300 MW and a2=0.1, rivals fixed at 50 MW each, Pi_1=0.2 q_1(200-q_1)+C_1. Moving q_1 from50 to100 MW gains $500/h. Joint adoption benefits do not prove unilateral stability.

**Response options / remaining work:** Required correction/clarification. Distinguish Nash equilibrium of prescribed J_i from monetary-profit maximization. Exact deviation calculation was added to earlier revised files, not the current original working source.

#### R2.3

**Review comment:** Variances and covariances enter the objectives only through terms independent of the bids. The same bidding rules therefore arise in a deterministic version of the model. This substantially limits the claimed contribution concerning uncertainty and risk-aware decisions.

**Issue:** The price-deviation metric decomposes into real-time price variance plus squared mean-price bias. At fixed means, coefficients, and bid limits, uncertainty contributes bid-independent constants to both expected profit and the objective. Hence optimal bids do not respond to covariance magnitude; the same rules arise in the deterministic case. The metric ignores individual imbalance exposure and loss direction and can be positive without uncertainty.

**Response options / remaining work:** Required clarification of interpretation and claim scope. Justify a quadratic settlement-price-deviation objective, not comprehensive downside risk protection. Changing to an uncertainty-responsive risk objective would change the model and may belong in future work. Sensitivity experiments can supplement the explanation on GitHub.

#### R2.4

**Review comment:** Equation (15) contains an incorrect covariance expansion.

**Issue:** Expected profit requires the identity E[x_i x_j]=E[x_i]E[x_j]+sigma_ij. Under the stated affine prices, its covariance contribution is 2a2(sigma_iL-sum_j sigma_ij). Incorrect signs or coefficients affect payoff levels.

**Response options / remaining work:** Required mathematical correction, not future work. Check downstream equations and numerical values. The correct contribution is bid-independent under this model, so explain why bid derivatives are unaffected. Independent moment calculations can be supplementary.

#### R2.5

**Review comment:** The dispatch equations impose no conventional generation bounds, and Gaussian generation and demand can produce physically infeasible realizations. The closed-form price expressions rely on this unconstrained setting, which is not realistic.

**Issue:** Unbounded bids and Gaussian physical quantities do not guarantee nonnegative, capacity-feasible conventional dispatch. Affine-price formulas assume an interior clearing regime; checking only equilibrium dispatch misses possible infeasible unilateral deviations.

**Response options / remaining work:** Required assumption and applicability clarification. Earlier revisions provide whole-box interior conditions, bounded moment-matched examples, and capacity checks at visited microgrid states. Binding conventional constraints and scarcity/congestion pricing need a different clearing analysis. Gaussian results must remain labeled unconstrained benchmarks.

### Reviewer 3

#### R3.1

**Review comment:** Since each operator can earn more by unilaterally deviating from the proposed strategy, the authors should explain why self-interested operators would follow it without enforcement.

**Issue:** Prescribed coefficients and consensus do not make ERS self-enforcing. Owners can earn more by overriding the controller. Merely assuming adherence leaves implementation unexplained.

**Response options / remaining work:** Unresolved modeling decision. Candidate budget-balanced transfers implement ERS for enrolled operators under enforceable settlement and fixed verified inputs. Free opt-out while others continue fails. Unanimous entry with a credible all-RN fallback can support participation, but requires binding commitments. This introduces an aggregator/coalition agreement and may strengthen the withholding objection. Do not insert automatically; institutional relevance must be justified.

#### R3.2

**Review comment:** The risk term depends on the bids only through a deterministic price spread, so the strategies are independent of uncertainty, and the authors should justify calling this parameter a risk attitude.

**Issue:** The price-deviation metric decomposes into real-time price variance plus squared mean-price bias. At fixed means, coefficients, and bid limits, uncertainty contributes bid-independent constants to both expected profit and the objective. Hence optimal bids do not respond to covariance magnitude; the same rules arise in the deterministic case. The metric ignores individual imbalance exposure and loss direction and can be positive without uncertainty.

**Response options / remaining work:** Required clarification of interpretation and claim scope. Justify a quadratic settlement-price-deviation objective, not comprehensive downside risk protection. Changing to an uncertainty-responsive risk objective would change the model and may belong in future work. Sensitivity experiments can supplement the explanation on GitHub.

#### R3.3

**Review comment:** The Pareto optimality is achieved by jointly withholding day-ahead bids to raise prices, so the authors should analyze its impact on consumers and market efficiency.

**Issue:** Maximizing renewable operators' aggregate profit is not social welfare maximization or Pareto improvement for all market participants. Lower day-ahead quantities can increase prices and consumer payments. Internal producer transfers do not fix this.

**Response options / remaining work:** Unresolved substantive analysis. Specify consumer settlement and conventional-generator accounting before calculating welfare effects. A prior audit derived a conditional baseline consumer-payment increase, but it is not an accepted manuscript result. Do not claim social efficiency from producer-profit results.

#### R3.4

**Review comment:** The distributed algorithm is flawed, since unequal consensus gains lead to a weighted rather than simple average, the utility directly receives each operator private forecast, and convergence slows as the number of operators grows, which the authors should correct.

**Issue:** Heterogeneous node gains generally fail to preserve the arithmetic mean. A connected graph alone is insufficient for the original averaging claim. Privacy and size-dependent convergence are additional separate problems.

**Response options / remaining work:** Required algorithmic correction. A common admissible gain on a connected undirected graph gives a symmetric doubly stochastic matrix. Check stopping criteria and estimation accuracy separately. Privacy needs its own observer assumptions and proof or a narrowed claim.

#### R3.5

**Review comment:** The case study only verifies the closed-form results on a small example, so the authors should add realistic studies with microgrid dynamics, network constraints, and comparisons with existing bidding methods.

**Issue:** Larger synthetic cases alone do not answer requests for network-constrained realism or comparisons with established bidding methods. Risk-neutral versus ERS is not a comparison with an external risk-sensitive benchmark.

**Response options / remaining work:** Unresolved evidence/scope decision. Select comparable existing methods and define fair information/objective constraints if adding a benchmark. Electrical network constraints require substantive modeling. These cannot be claimed resolved by export caps or aggregate balance alone.

### Reviewer 4

Recommendation: Major Revision.

#### R4.1

**Review comment:** The covariance terms in (15) appear inconsistent with (7a)–(7b) and Lemma 1. The authors should verify the signs and coefficients of these terms and clarify their impact on subsequent results.

**Issue:** Expected profit requires the identity E[x_i x_j]=E[x_i]E[x_j]+sigma_ij. Under the stated affine prices, its covariance contribution is 2a2(sigma_iL-sum_j sigma_ij). Incorrect signs or coefficients affect payoff levels.

**Response options / remaining work:** Required mathematical correction, not future work. Check downstream equations and numerical values. The correct contribution is bid-independent under this model, so explain why bid derivatives are unaffected. Independent moment calculations can be supplementary.

#### R4.2

**Review comment:** Under (7b), the real-time price variance is independent of the bidding coefficients. The authors should clarify the economic meaning of the proposed risk measure and its advantages over conventional bidding strategies.

**Issue:** The price-deviation metric decomposes into real-time price variance plus squared mean-price bias. At fixed means, coefficients, and bid limits, uncertainty contributes bid-independent constants to both expected profit and the objective. Hence optimal bids do not respond to covariance magnitude; the same rules arise in the deterministic case. The metric ignores individual imbalance exposure and loss direction and can be positive without uncertainty.

**Response options / remaining work:** Required clarification of interpretation and claim scope. Justify a quadratic settlement-price-deviation objective, not comprehensive downside risk protection. Changing to an uncertainty-responsive risk objective would change the model and may belong in future work. Sensitivity experiments can supplement the explanation on GitHub.

#### R4.3

**Review comment:** The proposed strategies may yield bidding coefficients exceeding unity, as shown in Table III. The authors should clarify the feasible bidding range and discuss the applicability of the theoretical results under practical operating constraints.

**Issue:** Uncapped equilibrium formulas can produce alpha_i>1. A finite bid constraint changes best responses and potentially the equilibrium and efficient coefficient design; clipping an uncapped solution is insufficient.

**Response options / remaining work:** Required feasible-set clarification and analysis if practical capped applicability is claimed. Earlier revisions contain constrained equilibrium and retuning results. Review for selective incorporation. Large random checks may be supplementary.

#### R4.4

**Review comment:** Lemma 2 permits heterogeneous consensus gains, which do not generally guarantee convergence to the arithmetic average. The authors should specify sufficient conditions for average consensus to ensure the correctness of Algorithm 1.

**Issue:** Heterogeneous node gains generally fail to preserve the arithmetic mean. A connected graph alone is insufficient for the original averaging claim. Privacy and size-dependent convergence are additional separate problems.

**Response options / remaining work:** Required algorithmic correction. A common admissible gain on a connected undirected graph gives a symmetric doubly stochastic matrix. Check stopping criteria and estimation accuracy separately. Privacy needs its own observer assumptions and proof or a narrowed claim.

#### R4.5

**Review comment:** Algorithm 1 initializes consensus variables with private generation forecasts, which may be disclosed through communication. The authors should clarify the privacy assumptions and justify the claimed privacy-preserving property.

**Issue:** Ordinary consensus is distributed but does not conceal forecasts from recipients. Algorithm 1 initializes psi_i(0)=x_i^o; a neighbor receiving that state receives the forecast. Computing a percentage locally at the end does not undo disclosure. On a star, a known-gain first updated leaf state can also reveal the forecast.

**Response options / remaining work:** Required claim correction or protocol change. Decide whom data must be hidden from. Either describe limited neighboring communication without a confidentiality guarantee, or add a protected aggregation protocol and proof. A GitHub replay can illustrate disclosure but cannot replace the argument.

#### R4.6

**Review comment:** The numerical studies are limited to a three-operator system, and explicit dynamic microgrid models are absent. The authors should strengthen the numerical validation and clarify the relevance to dynamic microgrid applications.

**Issue:** A static electricity-market game is not a model of electrical microgrid dynamics. Battery energy evolution demonstrates supervisory operation but does not establish voltage/frequency stability or network feasibility.

**Response options / remaining work:** Clarify application scope. Earlier synthetic sequential battery studies offer supporting evidence, but network power flow, fast dynamics, and strategic storage remain separate extensions unless explicitly added.

## Priority decisions

1. Privacy: narrow the claim or modify the protocol? Define who is trusted first.
2. Risk interpretation: retain price-deviation objective with clear limits, or change the model? No uncertainty adaptation is established at fixed inputs.
3. Incentives: retain a conditional controller-design result or justify an enforceable aggregator institution? An artificial transfer rule alone does not establish a realistic market arrangement.
4. Welfare: distinguish producer efficiency from consumer and total-market effects.
5. Scope and evidence: retain supervisory microgrid scope or add network/electrical modeling and external benchmarks?

## Supporting analysis already available locally

- `output/reviewer_notes/privacy_issue_review.md`: detailed first-comment discussion.
- `output/incentives/contracted_analysis.md`: candidate settlement proofs and entry assumptions.
- `output/incentives/optout_analysis.md`: explicit exit counterexamples and fallback distinctions.
- `output/incentives/skeptical_audit.md`: independent limitations, including coalition and consumer-accounting caveats.
- `output/dispatch/feasibility_checks.json`: prior interior-dispatch validation.
- Earlier numerical evidence is stored under `output/risk_caps`, `output/privacy`, `output/microgrid`, and `output/communication`.

These artifacts document prior investigation. They do not mean their assumptions or results have been approved for incorporation into the original manuscript.
