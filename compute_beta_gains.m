function K = compute_beta_gains(Sigma, beta, k_base, k_min, k_max, n_UAV)
% COMPUTE_BETA_GAINS  Individual risk-parameterized gain per agent.
%
%   k_i(beta) = k_base / (1 + beta * Sigma_i^2)
%
%   beta > 0  ->  gain reduced    (risk-averse:  cautious approach to orbit)
%   beta = 0  ->  gain = k_base  (risk-neutral:  standard consensus)
%   beta < 0  ->  gain amplified  (risk-seeking: aggressive orbit tracking)

denom = 1 + beta .* (Sigma .^ 2);
denom = max(denom, 0.05);
K     = k_base ./ denom;
K     = max(k_min, min(k_max, K));
end