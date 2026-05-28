%% Distributed Risk-Aware Nash Optimization Algorithm
% Implementation of Algorithm 1 from the paper:
% "Distributed Risk-Aware Bidding Strategy for Incorporating Renewable
%  Generation into Real-Time Electricity Market"
%
% ALGORITHM STRUCTURE (Algorithm 1 in paper):
%   Phase 1: Two consensus protocols (xi, psi) run until convergence
%            - Utility (node n+1): xi_{n+1}(0)=1, psi_{n+1}(0)=0
%            - Operator i:         xi_i(0)=0,      psi_i(0)=x_i^o
%   Phase 2: Each operator locally computes:
%            - n_hat = 1/xi_i - 1
%            - eta_i = x_i^o / (psi_i/xi_i)   [psi/xi -> x_r^o]
%            - beta_i' = (1-n_hat)/(4*n_hat)   [fair symmetric Pareto choice]
%   Phase 3: Submit bid:
%            - alpha_i* = 0.5*(1+4*beta_i') / eta_i    [eq. alpha_optimal_simplified]
%
% n is assumed known (per problem statement), but is also estimated via
% consensus for generality/verification.

clc; clear; close all;
set(groot, 'defaultFigureWindowStyle', 'docked');
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesColor', 'w');
set(groot, 'defaultAxesXColor', 'k');
set(groot, 'defaultAxesYColor', 'k');
set(groot, 'defaultTextColor', 'k');
set(groot, 'defaultLegendTextColor', 'k');
set(groot, 'defaultLegendColor', 'w');
set(groot, 'defaultLegendEdgeColor', 'k');
set(groot, 'defaultColorbarColor', 'k');

%% ========================================================================
%                           USER PARAMETERS
% =========================================================================

n  = 3;          % Number of renewable operators (known)
a2 = 0.1;        % Cost coefficients (utility)
a1 = 1;
a0 = 0;

x0  = [120; 100; 80];   % Forecasted renewable generation [MW]
xL0 = 400;              % Forecasted load [MW]

% Covariance (zero = deterministic)
sigma_r2 = 0;
sigma_L2 = 0;
sigma_rL = 0;

% Communication: fully connected graph over n operators + 1 utility
% Node ordering: [op1, op2, ..., opN, utility]
N_total = n + 1;
Adj = ones(N_total) - eye(N_total);   % fully connected
Deg = diag(sum(Adj, 2));
L_full = Deg - Adj;

% Consensus parameters
epsilon_c     = 0.1;      % step size (< 1/max_degree = 1/n)
max_cons_iter = 500;
cons_tol      = 1e-10;
plot_tail_buffer = 12;              % extra samples after detected convergence
cons_plot_tol = 1e-6;               % practical convergence tolerance for plotting horizon

% Algorithm 2 parameters (non-equal beta_i'(0))
beta0_non_equal = [-0.23; -0.08; 0.06];   % arbitrary initial beta_i'(0)
epsilon_beta    = 0.25;                    % outer update gain, require 0 < n*epsilon_beta < 2
delta_beta      = 1e-4;                    % keeps beta_i' > -1/4
beta_tol        = 1e-8;
max_beta_iter   = 250;
run_algo2       = false;                   % set true to re-enable Algorithm 2 section

% Algorithm 1 overbidding study (beta' sweep)
beta_sweep      = linspace(-0.24, 0.35, 220);

save_figs  = true;
plot_figs  = true;   % set false to skip all plotting (just run computations)
close_figs = true;   % set true to close all figures after saving
fig_out_dir = 'figures';

%% ========================================================================
%                         DERIVED QUANTITIES
% =========================================================================

xr0   = sum(x0);
eta   = x0 / xr0;           % true eta_i = x_i^o / x_r^o
fprintf('==========================================================\n');
fprintf('  DISTRIBUTED RISK-AWARE NASH OPTIMIZATION (Algorithm 1)\n');
fprintf('==========================================================\n\n');
fprintf('Parameters: n=%d, a2=%.2f, a1=%.2f, a0=%.2f\n', n, a2, a1, a0);
fprintf('x0 = [%s] MW,  xL0 = %.0f MW,  x_r^o = %.0f MW\n', num2str(x0'), xL0, xr0);
fprintf('eta = [%s]\n\n', num2str(eta', '%.6f '));

%% ========================================================================
%              PART 1: CENTRALIZED SOLUTION (GROUND TRUTH)
% =========================================================================

fprintf('----------------------------------------------------------\n');
fprintf('  PART 1: CENTRALIZED SOLUTION\n');
fprintf('----------------------------------------------------------\n\n');

% 1a. Coalition (Lemma 1): all operators cooperate, common alpha*=0.5, beta*=0
%     Achieves cooperative optimum g_max = 0.25*xr0^2  [eq. alpha_1 with beta=0]
alpha_coal = 0.5 * ones(n, 1);
f_coal     = sum(alpha_coal .* eta);        % = 0.5 since sum(eta)=1
g_coal     = f_coal * (1-f_coal) * xr0^2;  % = 0.25*xr0^2 = g_max
Pi_coal    = profit(alpha_coal, a2, a1, x0, xL0, xr0, eta, sigma_r2, sigma_L2, sigma_rL);
fprintf('Coalition / Lemma 1 (alpha*=0.5 for all, beta*=0):\n');
fprintf('  alpha = [%s]\n', num2str(alpha_coal', '%.6f '));
fprintf('  f = %.6f,  sum(Pi) = %.2f,  g = %.2f MW^2  (= g_max = %.2f MW^2)\n\n', ...
    f_coal, sum(Pi_coal), g_coal, 0.25*xr0^2);

% 1b. Risk-Neutral NE (beta_i' = 0)
[alpha_RN, f_RN, Pi_RN, ~] = nash_eq(zeros(n,1), n, a2, a1, x0, xL0, xr0, eta, sigma_r2, sigma_L2, sigma_rL);
fprintf('Risk-Neutral NE (beta''=0):\n');
fprintf('  alpha = [%s]\n', num2str(alpha_RN', '%.6f '));
fprintf('  f = %.6f,  sum(Pi) = %.2f,  g = %.2f MW^2  (g_max = %.2f MW^2)\n\n', ...
    f_RN, sum(Pi_RN), f_RN*(1-f_RN)*xr0^2, 0.25*xr0^2);

% 1c. Symmetric Pareto-Optimal: each beta_i' = (1-n)/(4n) => sum = (1-n)/4
beta_prime_sym = (1-n) / (4*n);
[alpha_sym, f_sym, Pi_sym, ~] = nash_eq(beta_prime_sym*ones(n,1), n, a2, a1, x0, xL0, xr0, eta, sigma_r2, sigma_L2, sigma_rL);
fprintf('Symmetric Pareto-Optimal [beta_i'' = (1-n)/(4n) = %.6f]:\n', beta_prime_sym);
fprintf('  alpha = [%s]\n', num2str(alpha_sym', '%.6f '));
fprintf('  f = %.6f,  sum(Pi) = %.2f,  g = %.2f MW^2  (g_max = %.2f MW^2)\n\n', ...
    f_sym, sum(Pi_sym), f_sym*(1-f_sym)*xr0^2, 0.25*xr0^2);

%% ========================================================================
%              PART 2: DISTRIBUTED ALGORITHM 1
% =========================================================================

fprintf('----------------------------------------------------------\n');
fprintf('  PART 2: DISTRIBUTED ALGORITHM 1 (One-Shot Consensus)\n');
fprintf('----------------------------------------------------------\n\n');

% ---- Phase 1: Consensus initialization ----
% xi:  utility=1, operators=0
% psi: utility=0, operators=x_i^o
xi_init  = zeros(N_total, 1);  xi_init(n+1)  = 1;
psi_init = zeros(N_total, 1);  psi_init(1:n) = x0;

fprintf('Phase 1: Running consensus protocols...\n');
fprintf('  xi_init  = [%s]\n', num2str(xi_init'));
fprintf('  psi_init = [%s]\n\n', num2str(psi_init'));

% Run both protocols and record full history for plotting
xi_hist  = zeros(N_total, max_cons_iter);
psi_hist = zeros(N_total, max_cons_iter);
xi  = xi_init;
psi = psi_init;
xi_iters  = max_cons_iter;
psi_iters = max_cons_iter;

for k = 1:max_cons_iter
    xi_hist(:,k)  = xi;
    psi_hist(:,k) = psi;
    xi_new  = xi  - epsilon_c * L_full * xi;
    psi_new = psi - epsilon_c * L_full * psi;
    if k > 1
        if norm(xi_new - xi) < cons_tol && xi_iters == max_cons_iter
            xi_iters = k;
        end
        if norm(psi_new - psi) < cons_tol && psi_iters == max_cons_iter
            psi_iters = k;
        end
    end
    xi  = xi_new;
    psi = psi_new;
end
xi_final  = xi;
psi_final = psi;

xi_bar  = mean(xi_final);
psi_bar = mean(psi_final);
fprintf('  xi  converged ~iter %d,  xi_bar  = %.8f  (true 1/%d = %.8f)\n', ...
    xi_iters, xi_bar, N_total, 1/N_total);
fprintf('  psi converged ~iter %d,  psi_bar = %.8f  (true x_r^o/%d = %.8f)\n\n', ...
    psi_iters, psi_bar, N_total, xr0/N_total);

% ---- Phase 2: Local computation at each operator ----
% Key insight: psi_i/xi_i -> (x_r^o/(n+1)) / (1/(n+1)) = x_r^o
% So each operator can recover x_r^o as the ratio, without knowing n+1 explicitly.
fprintf('Phase 2: Local computation at each operator:\n');

n_hat         = zeros(n, 1);
eta_dist      = zeros(n, 1);
xr0_dist      = zeros(n, 1);
beta_prime_d  = zeros(n, 1);

for i = 1:n
    xi_i  = xi_final(i);
    psi_i = psi_final(i);

    % Estimate n: Algorithm 1 line: n_hat = 1/xi_i - 1
    n_hat(i) = 1/xi_i - 1;

    % Estimate x_r^o: psi_i/xi_i -> x_r^o  (ratio of consensus values)
    xr0_dist(i) = psi_i / xi_i;

    % Estimate eta_i = x_i^o / x_r^o  (local x_i^o is known to operator i)
    eta_dist(i) = x0(i) / xr0_dist(i);

    % Fair symmetric Pareto beta: beta_i' = (1 - n_hat) / (4 * n_hat)
    % Ensures sum_i beta_i' = (1-n)/4  [eq. beta_sum_condition_simplified]
    beta_prime_d(i) = (1 - n_hat(i)) / (4 * n_hat(i));

    fprintf('  Op %d: n_hat=%.4f, x_r^o_est=%.4f, eta_est=%.6f, beta''=%.6f\n', ...
        i, n_hat(i), xr0_dist(i), eta_dist(i), beta_prime_d(i));
end

fprintf('\n  True:  n=%d, x_r^o=%.4f, eta=[%s]\n', n, xr0, num2str(eta', '%.6f '));
fprintf('  sum(beta'') = %.6f  (target (1-n)/4 = %.6f)\n\n', sum(beta_prime_d), (1-n)/4);

% ---- Phase 3: Bid submission ----
% alpha_i* = 0.5 * (1 + 4*beta_i') / eta_i   [eq. alpha_optimal_simplified]
fprintf('Phase 3: Computing and submitting bids...\n');
alpha_dist = 0.5 * (1 + 4*beta_prime_d) ./ eta_dist;
alpha_dist = max(0, min(1, alpha_dist));   % clamp to [0,1]

f_dist = sum(alpha_dist .* eta_dist);
g_dist = f_dist * (1-f_dist) * xr0^2;
Pi_dist = profit(alpha_dist, a2, a1, x0, xL0, xr0, eta_dist, sigma_r2, sigma_L2, sigma_rL);

for i = 1:n
    fprintf('  Op %d: alpha* = %.6f\n', i, alpha_dist(i));
end

%% ========================================================================
%              PART 3: OVERBIDDING SCENARIO STUDY (Corollary 5)
% =========================================================================
%
% Scenario 1 (current params): all eta_i >= 1/(2n), no overbidding at Pareto optimum.
% Scenario 2: x3_o = 30 MW => eta_3 = 30/250 = 0.12 < 1/6 = 1/(2n).
%             Op3 overbids at the symmetric Pareto beta' = (1-n)/(4n).
%             Corollary 5 gives the optimal constrained beta'* that recovers f=0.5.

fprintf('\n----------------------------------------------------------\n');
fprintf('  PART 3: OVERBIDDING SCENARIO STUDY\n');
fprintf('----------------------------------------------------------\n\n');

% Per-operator overbidding thresholds for Scenario 1
overbid_thresh = (2*eta_dist - 1) / 4;
fprintf('Scenario 1 (x3=80 MW, xr0=300):\n');
for i = 1:n
    fprintf('  Op%d: eta=%.4f  thresh=%.4f  (1/(2n)=%.4f) => %s\n', ...
        i, eta_dist(i), overbid_thresh(i), 1/(2*n), ...
        ternary(eta_dist(i) < 1/(2*n), 'OVERBIDS at Pareto', 'OK'));
end
fprintf('  Pareto beta''=%+.4f, Omega_s^P = empty, no overbidding.\n\n', beta_prime_sym);

% --- Overbidding scenario (heterogeneous beta, original system) ---
% beta = [-0.20, -0.15, -0.05]: Op3 (smallest eta) overbids under one-shot formula.
% sum(beta) = -0.40 → f_unc = 0.7 < 1, so profits remain positive.
% Comparison: unclamped overbidding vs no-overbidding (clamped).
beta_cV  = [-0.20; -0.15; -0.05];   % heterogeneous beta scenario
% Use same eta_dist from Algorithm 1 (estimated from consensus)
alpha_cV_unc = 0.5 * (1 + 4*beta_cV) ./ eta_dist;   % one-shot, unclamped
alpha_cV_clp = max(0, min(1, alpha_cV_unc));          % no-overbidding projection
Omega_cV = find(alpha_cV_unc > 1);

f_cV_unc  = sum(alpha_cV_unc .* eta_dist);
g_cV_unc  = f_cV_unc * (1 - f_cV_unc) * xr0^2;
Pi_cV_unc = profit(alpha_cV_unc, a2, a1, x0, xL0, xr0, eta_dist, sigma_r2, sigma_L2, sigma_rL);

f_cV_clp  = sum(alpha_cV_clp .* eta_dist);
g_cV_clp  = f_cV_clp * (1 - f_cV_clp) * xr0^2;
Pi_cV_clp = profit(alpha_cV_clp, a2, a1, x0, xL0, xr0, eta_dist, sigma_r2, sigma_L2, sigma_rL);

fprintf('Overbidding scenario (beta''=[%s]):\n', num2str(beta_cV', '%.2f '));
for i = 1:n
    fprintf('  Op%d: alpha_unc=%.4f  %s\n', i, alpha_cV_unc(i), ...
        ternary(alpha_cV_unc(i) > 1, '<= OVERBIDS', ''));
end
fprintf('  Omega_s = {Op%s}\n', num2str(Omega_cV'));
fprintf('  Unclamped: f=%.4f  g=%.1f MW^2  sum(Pi)=%.2f\n', f_cV_unc, g_cV_unc, sum(Pi_cV_unc));
fprintf('  Clamped:   f=%.4f  g=%.1f MW^2  sum(Pi)=%.2f  (g_max=%.1f)\n\n', ...
    f_cV_clp, g_cV_clp, sum(Pi_cV_clp), 0.25*xr0^2);

% --- Reference table ---
fprintf('Table values for LaTeX:\n');
fprintf('%-18s | %6s | %6s | %6s | %6s | %6s | %8s | %8s\n', ...
    'Scenario', 'al1', 'al2', 'al3', 'f', 'g(MW2)', 'sum(Pi)', 'note');
fprintf('%s\n', repmat('-',1,80));
fprintf('%-18s | %6.4f | %6.4f | %6.4f | %6.4f | %8.1f | %8.2f | %s\n', ...
    'S1 Pareto-NE', alpha_sym(1), alpha_sym(2), alpha_sym(3), f_sym, ...
    f_sym*(1-f_sym)*xr0^2, sum(Pi_sym), 'no overbid');
fprintf('%-18s | %6.4f | %6.4f | %6.4f | %6.4f | %8.1f | %8.2f | %s\n', ...
    'Unc. overbid', alpha_cV_unc(1), alpha_cV_unc(2), alpha_cV_unc(3), f_cV_unc, ...
    g_cV_unc, sum(Pi_cV_unc), sprintf('Op%s overbid', num2str(Omega_cV')));
fprintf('%-18s | %6.4f | %6.4f | %6.4f | %6.4f | %8.1f | %8.2f | %s\n', ...
    'No-overbid', alpha_cV_clp(1), alpha_cV_clp(2), alpha_cV_clp(3), f_cV_clp, ...
    g_cV_clp, sum(Pi_cV_clp), 'clamped');

%% ========================================================================
%              PART 3: DISTRIBUTED ALGORITHM 2 (NON-EQUAL beta_i'(0))
% =========================================================================

if run_algo2
fprintf('\n----------------------------------------------------------\n');
fprintf('  PART 4: DISTRIBUTED ALGORITHM 2 (Heterogeneous beta_i''(0))\n');
fprintf('----------------------------------------------------------\n\n');

if length(beta0_non_equal) ~= n
    error('beta0_non_equal must have length n.');
end

% Keep beta feasible with theorem bound beta_i' > -1/4
beta_k = max(beta0_non_equal(:), -0.25 + delta_beta);
B_target = (1 - n) / 4;

% Consensus over operators only for estimating mean(beta(k))
Adj_ops = ones(n) - eye(n);
L_ops   = diag(sum(Adj_ops, 2)) - Adj_ops;

beta_hist      = zeros(n, max_beta_iter);
Bsum_hist      = zeros(max_beta_iter, 1);
beta_err_hist  = zeros(max_beta_iter, 1);
alpha2_unc_hist = zeros(n, max_beta_iter);
f2_hist         = zeros(max_beta_iter, 1);
g2_hist         = zeros(max_beta_iter, 1);
sat2_hist       = zeros(max_beta_iter, 1);

fprintf('Initial beta''(0) = [%s]\n', num2str(beta_k', '%.6f '));
fprintf('Target sum(beta'') = %.6f\n\n', B_target);

beta_iters = max_beta_iter;
for k = 1:max_beta_iter
    beta_hist(:,k) = beta_k;

    % Inner consensus to estimate average(beta(k)); operators do not share full vector
    z = beta_k;
    for t = 1:max_cons_iter
        z_new = z - epsilon_c * L_ops * z;
        if norm(z_new - z) < cons_tol
            z = z_new;
            break;
        end
        z = z_new;
    end

    B_est = n * mean(z);
    err_B = B_est - B_target;
    Bsum_hist(k) = B_est;
    beta_err_hist(k) = err_B;

    % Track one-shot bids induced by current beta_i'(k)
    alpha_unc_k = 0.5 * (1 + 4*beta_k) ./ eta_dist;
    alpha2_unc_hist(:,k) = alpha_unc_k;
    f2_hist(k) = sum(alpha_unc_k .* eta_dist);
    g2_hist(k) = f2_hist(k) * (1 - f2_hist(k)) * xr0^2;
    sat2_hist(k) = sum(alpha_unc_k > 1);

    if abs(err_B) < beta_tol
        beta_iters = k;
        break;
    end

    % Common correction preserves heterogeneity while steering sum(beta')
    beta_k = max(beta_k - epsilon_beta * err_B, -0.25 + delta_beta);
end

beta_hist = beta_hist(:,1:beta_iters);
Bsum_hist = Bsum_hist(1:beta_iters);
beta_err_hist = beta_err_hist(1:beta_iters);
alpha2_unc_hist = alpha2_unc_hist(:,1:beta_iters);
f2_hist = f2_hist(1:beta_iters);
g2_hist = g2_hist(1:beta_iters);
sat2_hist = sat2_hist(1:beta_iters);
beta_prime_algo2 = beta_k;

% One-shot bidding after sum(beta') coordination (unclamped assumption)
alpha_algo2 = 0.5 * (1 + 4*beta_prime_algo2) ./ eta_dist;
f_algo2 = sum(alpha_algo2 .* eta_dist);
g_algo2 = f_algo2 * (1-f_algo2) * xr0^2;
Pi_algo2 = profit(alpha_algo2, a2, a1, x0, xL0, xr0, eta_dist, sigma_r2, sigma_L2, sigma_rL);

fprintf('Converged in %d outer iterations.\n', beta_iters);
fprintf('Final beta'' = [%s]\n', num2str(beta_prime_algo2', '%.6f '));
fprintf('sum(beta'') = %.8f  (target %.8f)\n', sum(beta_prime_algo2), B_target);
fprintf('Saturated operators (alpha_unc > 1) at final iterate: %d/%d\n', sat2_hist(end), n);
for i = 1:n
    fprintf('  Op %d: alpha* = %.6f\n', i, alpha_algo2(i));
end
fprintf('  f = %.6f,  g = %.2f,  sum(Pi) = %.2f\n', f_algo2, g_algo2, sum(Pi_algo2));
end

%% ========================================================================
%         PART 4: CONSTRAINED RISK-NEUTRAL NASH (Corollary 4)
%
% With modified capacities (x3=30 MW, xr0=250), eta_3=0.12 < 1/(n+1)=0.25.
% Operator 3 would overbid under unconstrained RN. Corollary 4 gives the
% constrained RN solution by clamping Op3 to alpha_3=1 and re-solving.
% Compare against: (a) unconstrained RN, (b) Cor-4 clamped RN,
%                  (c) ERS Pareto via heterogeneous-beta partition (Cor. 5b).
% =========================================================================

fprintf('\n----------------------------------------------------------\n');
fprintf('  PART 4: CONSTRAINED RN NASH — COROLLARY 4\n');
fprintf('----------------------------------------------------------\n\n');

x0_s2   = [120; 100; 30];          % Scenario 2 capacities
xr0_s2  = sum(x0_s2);              % 250 MW
eta_s2  = x0_s2 / xr0_s2;         % [0.48, 0.40, 0.12]

fprintf('Scenario 2 params: x0=[%s] MW,  xr0=%.0f MW\n', num2str(x0_s2'), xr0_s2);
fprintf('eta=[%s],  RN threshold 1/(n+1)=%.4f\n\n', num2str(eta_s2', '%.4f '), 1/(n+1));

% (a) Unconstrained RN (beta=0, overbidding allowed)
[alpha_RN_s2, f_RN_s2, Pi_RN_s2, ~] = nash_eq(zeros(n,1), n, a2, a1, x0_s2, xL0, xr0_s2, eta_s2, 0,0,0);
g_RN_s2 = f_RN_s2 * (1-f_RN_s2) * xr0_s2^2;
fprintf('(a) Unconstrained RN (beta=0):\n');
fprintf('    alpha=[%s]\n', num2str(alpha_RN_s2', '%.4f '));
fprintf('    f=%.4f  g=%.1f MW^2  sum(Pi)=%.2f  (g_max=%.0f)\n\n', ...
    f_RN_s2, g_RN_s2, sum(Pi_RN_s2), 0.25*xr0_s2^2);

% (b) Corollary 4: constrained RN (beta=0, no overbidding)
Omega_RN_s2 = find(eta_s2 < 1/(n+1));          % operators forced to clamp
S_RN_s2     = sum(eta_s2(Omega_RN_s2));         % sum of clamped etas
n_eff_RN    = n + 1 - length(Omega_RN_s2);      % effective n+1-|Omega|
alpha_Cor4  = zeros(n, 1);
for i = 1:n
    if ismember(i, Omega_RN_s2)
        alpha_Cor4(i) = 1;
    else
        alpha_Cor4(i) = (1/n_eff_RN) * (1/eta_s2(i)) * (1 - S_RN_s2);
    end
end
f_Cor4  = sum(alpha_Cor4 .* eta_s2);
g_Cor4  = f_Cor4 * (1-f_Cor4) * xr0_s2^2;
Pi_Cor4 = profit(alpha_Cor4, a2, a1, x0_s2, xL0, xr0_s2, eta_s2, 0,0,0);
fprintf('(b) Corollary 4 — constrained RN (Omega_RN={Op%s}):\n', num2str(Omega_RN_s2'));
fprintf('    alpha=[%s]\n', num2str(alpha_Cor4', '%.4f '));
fprintf('    f=%.4f  g=%.1f MW^2  sum(Pi)=%.2f\n\n', f_Cor4, g_Cor4, sum(Pi_Cor4));

% (c) Corollary 5b — Heterogeneous-beta partition (Scenario 2)
%   eta_c = (eta_max + eta_min)/2; Omega={i: eta_i>=eta_c}; Omega^c={i: eta_i<eta_c}
%   For i in Omega^c: alpha_i*=1.  For i in Omega: alpha_i*=D/(2*|Omega|*eta_i)
%   NOTE: For S2 this gives numerically identical result to old "Remark 1"
%         (the partitions coincide: Omega^c = Omega_ERS).
eta_max_s2   = max(eta_s2);
eta_min_s2   = min(eta_s2);
eta_c_s2     = (eta_max_s2 + eta_min_s2) / 2;
Omega_ERS_s2  = find(eta_s2 >= eta_c_s2);     % large operators (in Omega)
OmegaC_ERS_s2 = find(eta_s2 <  eta_c_s2);    % small operators (in Omega^c, alpha=1)
S_ERS_s2     = sum(eta_s2(Omega_ERS_s2));
S_ERS_C_s2   = sum(eta_s2(OmegaC_ERS_s2));
D_s2         = S_ERS_s2 - S_ERS_C_s2;
nOmega_s2    = length(Omega_ERS_s2);
alpha_Rem1   = zeros(n, 1);     % kept as "alpha_Rem1" for downstream compatibility
for i = 1:n
    if ismember(i, OmegaC_ERS_s2)
        alpha_Rem1(i) = 1;
    else
        alpha_Rem1(i) = D_s2 / (2 * nOmega_s2 * eta_s2(i));
    end
end
f_Rem1  = sum(alpha_Rem1 .* eta_s2);
g_Rem1  = f_Rem1 * (1-f_Rem1) * xr0_s2^2;
Pi_Rem1 = profit(alpha_Rem1, a2, a1, x0_s2, xL0, xr0_s2, eta_s2, 0,0,0);
fprintf('(c) Corollary 5b — Heterogeneous-beta partition (eta_c=%.3f, Omega={Op%s}, Omega^c={Op%s}):\n', ...
    eta_c_s2, num2str(Omega_ERS_s2'), num2str(OmegaC_ERS_s2'));
fprintf('    D=%.4f  |Omega|=%d\n', D_s2, nOmega_s2);
fprintf('    alpha=[%s]\n', num2str(alpha_Rem1', '%.4f '));
fprintf('    f=%.4f  g=%.1f MW^2  sum(Pi)=%.2f  (g_max=%.0f)\n\n', ...
    f_Rem1, g_Rem1, sum(Pi_Rem1), 0.25*xr0_s2^2);

%% ========================================================================
%         PART 5: HETEROGENEOUS PARETO ALLOCATIONS (Theorem 2)
%
% Theorem 2 states that *any* beta with sum(beta') = (1-n)/4 achieves Pareto
% optimality. Here we demonstrate three distinct allocations for n=3
% (original capacities) that all yield f=0.5, g=g_max, but distribute
% individual profits differently. This illustrates the degrees of freedom
% available to a regulator or cooperative agreement.
%   Target: sum(beta') = (1-3)/4 = -0.5
% =========================================================================

fprintf('----------------------------------------------------------\n');
fprintf('  PART 5: HETEROGENEOUS PARETO ALLOCATIONS — THEOREM 2\n');
fprintf('----------------------------------------------------------\n\n');

% Three valid heterogeneous risk allocations (all sum to -0.5)
beta_het = [ ...
    -1/6,  -1/6,  -1/6;   % Row 1: symmetric ERS (baseline)
    -0.10, -0.17, -0.23;  % Row 2: Op1 least risk-seeking (larger share)
    -0.23, -0.17, -0.10;  % Row 3: Op1 most risk-seeking
];
het_labels = {'Symmetric ERS', 'Hetero A (\beta''=[-.10,-.17,-.23])', ...
              'Hetero B (\beta''=[-.23,-.17,-.10])'};

n_het = size(beta_het, 1);
alpha_het = zeros(n, n_het);
f_het     = zeros(1, n_het);
g_het     = zeros(1, n_het);
Pi_het    = zeros(n, n_het);

for k = 1:n_het
    bk = beta_het(k, :)';
    % Theorem 2 closed-form: alpha_i* = 0.5*(1+4*beta_i')/eta_i
    alpha_k = 0.5 * (1 + 4*bk) ./ eta;
    % Clamp to [0,1] (should not bind for these allocations with original x0)
    alpha_k = max(0, min(1, alpha_k));
    f_k     = sum(alpha_k .* eta);
    g_k     = f_k * (1-f_k) * xr0^2;
    Pi_k    = profit(alpha_k, a2, a1, x0, xL0, xr0, eta, 0,0,0);
    alpha_het(:,k) = alpha_k;
    f_het(k)       = f_k;
    g_het(k)       = g_k;
    Pi_het(:,k)    = Pi_k;
    fprintf('Allocation %d — %s:\n', k, het_labels{k});
    fprintf('  beta''=[%s]  sum=%.4f\n', num2str(bk', '%.4f '), sum(bk));
    fprintf('  alpha=[%s]\n', num2str(alpha_k', '%.4f '));
    fprintf('  f=%.4f  g=%.1f MW^2  sum(Pi)=%.2f\n', f_k, g_k, sum(Pi_k));
    fprintf('  Per-op Pi=[%s]\n\n', num2str(Pi_k', '%.2f '));
end
fprintf('Key insight: all three achieve f=0.5, g=g_max=%.0f MW^2,\n', 0.25*xr0^2);
fprintf('but individual profits differ — regulator flexibility.\n\n');

%% ========================================================================
%         PART 6: MULTIPLE SATURATION — COROLLARY 5b (|Omega_ERS|=2)
%
% To test generality of Corollary 5b beyond the single-operator case,
% we set x1=200, x2=50, x3=50 MW. Both Op2 and Op3 have
% eta_2=eta_3=50/300=0.167 = 1/(2n), sitting at the overbid boundary.
% We use slightly smaller values to ensure saturation: x2=x3=45 => eta=0.15.
% =========================================================================

fprintf('----------------------------------------------------------\n');
fprintf('  PART 6: MULTIPLE SATURATION — COROLLARY 5b (HETEROGENEOUS-BETA PARTITION)\n');
fprintf('----------------------------------------------------------\n\n');

x0_s6   = [210; 45; 45];
xr0_s6  = sum(x0_s6);              % 300 MW (same xr0 for comparability)
eta_s6  = x0_s6 / xr0_s6;         % [0.70, 0.15, 0.15]
thresh_s6 = 1/(2*n);               % 0.1667

fprintf('Params: x0=[%s] MW,  xr0=%.0f MW\n', num2str(x0_s6'), xr0_s6);
fprintf('eta=[%s],  ERS threshold 1/(2n)=%.4f\n', num2str(eta_s6', '%.4f '), thresh_s6);
fprintf('Omega_ERS = {');
Omega_s6 = find(eta_s6 < thresh_s6);
fprintf('Op%d ', Omega_s6); fprintf('}  (|Omega|=%d)\n\n', length(Omega_s6));

S_s6     = sum(eta_s6(Omega_s6));
n_eff_s6 = 2*n - length(Omega_s6);

% Unconstrained ERS (symmetric Pareto beta, no clamp)
beta_sym_s6  = (1-n)/(4*n);
alpha_unc_s6 = 0.5 * (1 + 4*beta_sym_s6) ./ eta_s6;
f_unc_s6     = sum(alpha_unc_s6 .* eta_s6);
g_unc_s6     = f_unc_s6 * (1-f_unc_s6) * xr0_s6^2;
Pi_unc_s6    = profit(alpha_unc_s6, a2, a1, x0_s6, xL0, xr0_s6, eta_s6, 0,0,0);
fprintf('Unconstrained ERS (beta''=%.4f, allows overbidding):\n', beta_sym_s6);
fprintf('  alpha=[%s]\n', num2str(alpha_unc_s6', '%.4f '));
fprintf('  f=%.4f  g=%.1f MW^2  sum(Pi)=%.2f\n\n', f_unc_s6, g_unc_s6, sum(Pi_unc_s6));

% Corollary 5b — Heterogeneous-beta partition construction
%   eta_c = (eta_max + eta_min)/2
%   Omega  = {i : eta_i >= eta_c}   (large-capacity operators)
%   Omega^c = {i : eta_i <  eta_c}  (small-capacity, alpha_i*=1)
%   For i in Omega: alpha_i* = D / (2*|Omega|*eta_i), D = S_Omega - S_OmegaC
eta_max_s6   = max(eta_s6);
eta_min_s6   = min(eta_s6);
eta_c_s6     = (eta_max_s6 + eta_min_s6) / 2;
Omega_het_s6  = find(eta_s6 >= eta_c_s6);   % {1}
OmegaC_het_s6 = find(eta_s6 <  eta_c_s6);  % {2, 3}
S_Omega_s6    = sum(eta_s6(Omega_het_s6));
S_OmegaC_s6  = sum(eta_s6(OmegaC_het_s6));
D_s6         = S_Omega_s6 - S_OmegaC_s6;
nOmega_s6    = length(Omega_het_s6);

alpha_het_s6 = zeros(n, 1);
for i = 1:n
    if ismember(i, OmegaC_het_s6)
        alpha_het_s6(i) = 1;
    else
        alpha_het_s6(i) = D_s6 / (2 * nOmega_s6 * eta_s6(i));
    end
end
f_het_s6  = sum(alpha_het_s6 .* eta_s6);
g_het_s6  = f_het_s6 * (1-f_het_s6) * xr0_s6^2;
Pi_het_s6 = profit(alpha_het_s6, a2, a1, x0_s6, xL0, xr0_s6, eta_s6, 0,0,0);
fprintf('Corollary 5b — Heterogeneous-beta partition (eta_c=%.3f, Omega={Op%s}):\n', ...
    eta_c_s6, num2str(Omega_het_s6'));
fprintf('  D = S_Omega-S_OmegaC = %.4f-%.4f = %.4f,  |Omega|=%d\n', ...
    S_Omega_s6, S_OmegaC_s6, D_s6, nOmega_s6);
fprintf('  alpha=[%s]\n', num2str(alpha_het_s6', '%.4f '));
fprintf('  f=%.4f  g=%.1f MW^2  sum(Pi)=%.2f  (g_max=%.0f)\n\n', ...
    f_het_s6, g_het_s6, sum(Pi_het_s6), 0.25*xr0_s6^2);

%% ========================================================================
%    PART 7: DISTRIBUTED IMPLEMENTATIONS OF NEW SCENARIOS
%
%  7a: Hetero A — Algorithm 1 consensus + regulator-assigned beta_i'
%  7b: Constrained ERS Scenario 2 (|Omega|=1) — third phi consensus
%  7c: Multiple saturation Scenario 3 (|Omega|=2) — third phi consensus
% =========================================================================

fprintf('\n----------------------------------------------------------\n');
fprintf('  PART 7: DISTRIBUTED IMPLEMENTATIONS\n');
fprintf('----------------------------------------------------------\n\n');

% ---- 7a: Distributed Heterogeneous Pareto (Hetero A) --------------------
% Algorithm 1 Phase 1+2 already ran (eta_dist, n_hat from base x0).
% Each operator uses a regulator-assigned beta_i' rather than the symmetric
% formula in Phase 2; everything else in Algorithm 1 is unchanged.
fprintf('7a — Distributed Hetero A  (beta''=[%.2f,%.2f,%.2f]):\n', ...
    beta_het(2,1), beta_het(2,2), beta_het(2,3));
beta_hetA_dist   = beta_het(2,:)';               % regulator-assigned
alpha_hetA_dist  = 0.5*(1+4*beta_hetA_dist)./eta_dist;
alpha_hetA_dist  = max(0, min(1, alpha_hetA_dist));
f_hetA_dist      = sum(alpha_hetA_dist.*eta_dist);
g_hetA_dist      = f_hetA_dist*(1-f_hetA_dist)*xr0^2;
Pi_hetA_dist     = profit(alpha_hetA_dist,a2,a1,x0,xL0,xr0,eta_dist,0,0,0);
fprintf('  alpha=[%s]  f=%.4f  g=%.1f  sum(Pi)=%.2f\n', ...
    num2str(alpha_hetA_dist','%.4f '), f_hetA_dist, g_hetA_dist, sum(Pi_hetA_dist));
fprintf('  Centralized error: %.2e\n\n', norm(alpha_hetA_dist - alpha_het(:,2)));

% ---- 7b: Distributed Constrained ERS — Scenario 2 (Omega^c={3}) ----------
fprintf('7b — Distributed Cor5b (hetero-beta), Scenario 2  (x3=30 MW):\n');

% Phase 1: xi / psi consensus with x0_s2
xi_init_s2  = zeros(N_total,1);  xi_init_s2(n+1)  = 1;
psi_init_s2 = zeros(N_total,1);  psi_init_s2(1:n) = x0_s2;
[xi_f_s2,  xi_hist_s2,  ~] = run_cons(xi_init_s2,  L_full, epsilon_c, max_cons_iter, cons_tol);
[psi_f_s2, psi_hist_s2, ~] = run_cons(psi_init_s2, L_full, epsilon_c, max_cons_iter, cons_tol);

% Phase 2: each operator locally estimates n and eta
n_hat_s2d = zeros(n,1);  eta_d_s2 = zeros(n,1);
for i = 1:n
    n_hat_s2d(i) = 1/xi_f_s2(i) - 1;
    eta_d_s2(i)  = x0_s2(i) / (psi_f_s2(i)/xi_f_s2(i));
end
n_est_s2 = round(mean(n_hat_s2d));

% Phase 3a: max-/min-consensus on eta (exact in 2 hops on star K_{1,n})
%   Each operator knows its own eta; utility collects all and broadcasts
%   max/min back.  Simulated here by global max/min of converged estimates.
eta_max_d_s2 = max(eta_d_s2);
eta_min_d_s2 = min(eta_d_s2);
eta_c_d_s2   = (eta_max_d_s2 + eta_min_d_s2) / 2;
% in_OmegaERS_d_s2: true = in Omega_ERS (small, eta < eta_c), gets alpha=1
% not_in_OmegaERS_d_s2: true = NOT in Omega_ERS (large, eta >= eta_c), gets formula
in_OmegaERS_d_s2     = eta_d_s2 <  eta_c_d_s2;
not_in_OmegaERS_d_s2 = ~in_OmegaERS_d_s2;

% Phase 4: xi2-consensus to recover (n - |Omega_ERS|)
%   Seed: xi2_i(0) = 1 if NOT in Omega_ERS, 0 if in Omega_ERS; utility seeds 0
%   After convergence: (n - |Omega_ERS|) = xi2_f / xi_f  (ratio trick)
xi2_init_s2 = zeros(N_total,1);
xi2_init_s2(1:n) = double(not_in_OmegaERS_d_s2);
[xi2_f_s2, ~, ~] = run_cons(xi2_init_s2, L_full, epsilon_c, max_cons_iter, cons_tol);
n_minus_nERS_s2 = round(xi2_f_s2(1) / xi_f_s2(1));

% Phase 4: psi2-consensus to recover sum_{j not in Omega_ERS} eta_j
%   Seed: psi2_i(0) = eta_i if NOT in Omega_ERS, 0 if in Omega_ERS; utility seeds 0
%   After convergence: sum_notERS = psi2_f / xi_f  (ratio trick)
psi2_init_s2 = zeros(N_total,1);
psi2_init_s2(1:n) = eta_d_s2 .* double(not_in_OmegaERS_d_s2);
[psi2_f_s2, ~, ~] = run_cons(psi2_init_s2, L_full, epsilon_c, max_cons_iter, cons_tol);
sum_notERS_s2 = psi2_f_s2(1) / xi_f_s2(1);

% Corollary 5b bids (distributed, Phase 4)
alpha_opt_s2d = zeros(n,1);
for i = 1:n
    if in_OmegaERS_d_s2(i)
        alpha_opt_s2d(i) = 1;
    else
        alpha_opt_s2d(i) = (2*sum_notERS_s2 - 1) / (2 * n_minus_nERS_s2 * eta_d_s2(i));
    end
end

fprintf('  eta_c=%.4f  not-Omega_ERS={Op%s}  Omega_ERS={Op%s}\n', ...
    eta_c_d_s2, num2str(find(not_in_OmegaERS_d_s2)'), num2str(find(in_OmegaERS_d_s2)'));
fprintf('  sum_notERS_hat=%.4f (true=%.4f)  (n-|OmegaERS|)_hat=%d (true=%d)\n', ...
    sum_notERS_s2, S_ERS_s2, n_minus_nERS_s2, nOmega_s2);
fprintf('  Cor5b dist: alpha=[%s]  err vs central=%.2e\n\n', ...
    num2str(alpha_opt_s2d','%.4f '), norm(alpha_opt_s2d - alpha_Rem1));

% ---- 7c: Distributed Constrained ERS — Scenario 3 (Omega^c={2,3}) -------
fprintf('7c — Distributed Cor5b (hetero-beta), Scenario 3  (x=[210,45,45] MW):\n');

% Phase 1
xi_init_s6  = zeros(N_total,1);  xi_init_s6(n+1)  = 1;
psi_init_s6 = zeros(N_total,1);  psi_init_s6(1:n) = x0_s6;
[xi_f_s6,  xi_hist_s6,  ~] = run_cons(xi_init_s6,  L_full, epsilon_c, max_cons_iter, cons_tol);
[psi_f_s6, psi_hist_s6, ~] = run_cons(psi_init_s6, L_full, epsilon_c, max_cons_iter, cons_tol);

% Phase 2: estimate n and eta
n_hat_s6d = zeros(n,1);  eta_d_s6 = zeros(n,1);
for i = 1:n
    n_hat_s6d(i) = 1/xi_f_s6(i) - 1;
    eta_d_s6(i)  = x0_s6(i) / (psi_f_s6(i)/xi_f_s6(i));
end
n_est_s6 = round(mean(n_hat_s6d));

% Phase 3a: max-/min-consensus on eta (exact in 2 hops on star K_{1,n})
eta_max_d_s6 = max(eta_d_s6);
eta_min_d_s6 = min(eta_d_s6);
eta_c_d_s6   = (eta_max_d_s6 + eta_min_d_s6) / 2;
% in_OmegaERS_d_s6: true = in Omega_ERS (small, eta < eta_c), gets alpha=1
% not_in_OmegaERS_d_s6: true = NOT in Omega_ERS (large, eta >= eta_c), gets formula
in_OmegaERS_d_s6     = eta_d_s6 <  eta_c_d_s6;
not_in_OmegaERS_d_s6 = ~in_OmegaERS_d_s6;

% Phase 4: xi2-consensus to recover (n - |Omega_ERS|)
%   Seed: xi2_i(0) = 1 if NOT in Omega_ERS, 0 if in Omega_ERS; utility seeds 0
%   After convergence: (n - |Omega_ERS|) = xi2_f / xi_f  (ratio trick)
xi2_init_s6 = zeros(N_total,1);
xi2_init_s6(1:n) = double(not_in_OmegaERS_d_s6);
[xi2_f_s6, ~, ~] = run_cons(xi2_init_s6, L_full, epsilon_c, max_cons_iter, cons_tol);
n_minus_nERS_s6 = round(xi2_f_s6(1) / xi_f_s6(1));

% Phase 4: psi2-consensus to recover sum_{j not in Omega_ERS} eta_j
%   Seed: psi2_i(0) = eta_i if NOT in Omega_ERS, 0 if in Omega_ERS; utility seeds 0
%   After convergence: sum_notERS = psi2_f / xi_f  (ratio trick)
psi2_init_s6 = zeros(N_total,1);
psi2_init_s6(1:n) = eta_d_s6 .* double(not_in_OmegaERS_d_s6);
[psi2_f_s6, ~, ~] = run_cons(psi2_init_s6, L_full, epsilon_c, max_cons_iter, cons_tol);
sum_notERS_s6 = psi2_f_s6(1) / xi_f_s6(1);

% Corollary 5b bids (distributed, Phase 4)
alpha_opt_s6d = zeros(n,1);
for i = 1:n
    if in_OmegaERS_d_s6(i)
        alpha_opt_s6d(i) = 1;
    else
        alpha_opt_s6d(i) = (2*sum_notERS_s6 - 1) / (2 * n_minus_nERS_s6 * eta_d_s6(i));
    end
end

fprintf('  eta_c=%.4f  not-Omega_ERS={Op%s}  Omega_ERS={Op%s}\n', ...
    eta_c_d_s6, num2str(find(not_in_OmegaERS_d_s6)'), num2str(find(in_OmegaERS_d_s6)'));
fprintf('  sum_notERS_hat=%.4f (true=%.4f)  (n-|OmegaERS|)_hat=%d (true=%d)\n', ...
    sum_notERS_s6, S_Omega_s6, n_minus_nERS_s6, nOmega_s6);
fprintf('  Cor5b dist: alpha=[%s]  err vs central=%.2e\n\n', ...
    num2str(alpha_opt_s6d','%.4f '), norm(alpha_opt_s6d - alpha_het_s6));

%% ========================================================================
%                         FINAL SUMMARY
% =========================================================================

fprintf('\n==========================================================\n');
fprintf('  FINAL RESULTS SUMMARY\n');
fprintf('==========================================================\n');
fprintf('%-35s %-14s %-14s %-14s %-14s\n', '', 'Coalition', 'Risk-Neutral', 'Pareto(Central)', 'Algo1(Dist)');
fprintf('%s\n', repmat('-', 1, 100));
for i = 1:n
    fprintf('  alpha_%d*                           %-14.6f %-14.6f %-14.6f %-14.6f\n', ...
        i, alpha_coal(i), alpha_RN(i), alpha_sym(i), alpha_dist(i));
end
fprintf('  f                                  %-14.6f %-14.6f %-14.6f %-14.6f\n', f_coal, f_RN, f_sym, f_dist);
fprintf('  g = f(1-f)*x_r^o^2 [MW^2]         %-14.2f %-14.2f %-14.2f %-14.2f\n', ...
    g_coal, f_RN*(1-f_RN)*xr0^2, f_sym*(1-f_sym)*xr0^2, g_dist);
fprintf('  sum(Pi)                            %-14.2f %-14.2f %-14.2f %-14.2f\n', ...
    sum(Pi_coal), sum(Pi_RN), sum(Pi_sym), sum(Pi_dist));
for i = 1:n
    fprintf('  Pi_%d                               %-14.2f %-14.2f %-14.2f %-14.2f\n', ...
        i, Pi_coal(i), Pi_RN(i), Pi_sym(i), Pi_dist(i));
end
fprintf('  Improvement vs Risk-Neutral        %-14.2f%% %-14s %-14.2f%% %-14.2f%%\n', ...
    100*(sum(Pi_coal)-sum(Pi_RN))/abs(sum(Pi_RN)), '-', ...
    100*(sum(Pi_sym)-sum(Pi_RN))/abs(sum(Pi_RN)), ...
    100*(sum(Pi_dist)-sum(Pi_RN))/abs(sum(Pi_RN)));
fprintf('  Pareto efficiency (g/g_max)        %-14.2f%% %-14.2f%% %-14.2f%% %-14.2f%%\n', ...
    100*g_coal/(0.25*xr0^2), 100*f_RN*(1-f_RN)*xr0^2/(0.25*xr0^2), ...
    100*f_sym*(1-f_sym)*xr0^2/(0.25*xr0^2), 100*g_dist/(0.25*xr0^2));
fprintf('\n  Overbidding scenario (beta''=[%.2f,%.2f,%.2f]):\n', beta_cV(1), beta_cV(2), beta_cV(3));
fprintf('    unclamped sum(Pi)               = %.2f  (f=%.3f, g=%.1f MW^2)\n', sum(Pi_cV_unc), f_cV_unc, g_cV_unc);
fprintf('    no-overbidding sum(Pi)          = %.2f  (f=%.3f, g=%.1f MW^2)\n', sum(Pi_cV_clp), f_cV_clp, g_cV_clp);

%% ========================================================================
%                         VISUALIZATION
% =========================================================================

if ~plot_figs
    fprintf('\nplot_figs=false: skipping all figures.\n');
    return;
end

op_labels = arrayfun(@(i) sprintf('Op%d', i), 1:n, 'UniformOutput', false);
node_labels = [op_labels, {'Utility'}];

% Practical plot horizon: first index where trajectories are close to final consensus.
xi_err = max(abs(xi_hist - xi_final), [], 1);
psi_err = max(abs(psi_hist - psi_final), [], 1);
k_xi_plot = find(xi_err < cons_plot_tol, 1, 'first');
k_psi_plot = find(psi_err < cons_plot_tol, 1, 'first');
if isempty(k_xi_plot), k_xi_plot = xi_iters; end
if isempty(k_psi_plot), k_psi_plot = psi_iters; end

K_plot = min(max_cons_iter, max([k_xi_plot, k_psi_plot]) + plot_tail_buffer);
K = max(2, K_plot);
% Plot A1: xi consensus trajectories
figure('Position', [50, 50, 620, 460]);
plot(1:K, xi_hist(1:n,1:K)', 'LineWidth', 1.5); hold on;
plot(1:K, xi_hist(n+1,1:K)', 'k--', 'LineWidth', 1.5);
yline(1/N_total, 'r:', 'LineWidth', 1.5);
xlabel('Iteration k'); ylabel('\xi_i(k)');
title('Consensus \xi (counting)');
legend([node_labels, {sprintf('1/(n+1)=%.4f',1/N_total)}], 'Location','best','FontSize',7);
grid on;
if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end
    fig1a_path = fullfile(fig_out_dir, 'algorithm1_xi_consensus.png');
    exportgraphics(gcf, fig1a_path, 'Resolution', 150);
    fprintf('\nFigure saved to: %s\n', fig1a_path);
end

% Plot A2: psi consensus trajectories
figure('Position', [80, 80, 620, 460]);
plot(1:K, psi_hist(1:n,1:K)', 'LineWidth', 1.5); hold on;
plot(1:K, psi_hist(n+1,1:K)', 'k--', 'LineWidth', 1.5);
yline(xr0/N_total, 'r:', 'LineWidth', 1.5);
xlabel('Iteration k'); ylabel('\psi_i(k)  (MW)');
title('Consensus \psi (forecast diffusion)');
legend([node_labels, {sprintf('x_r^o/(n+1)=%.1f',xr0/N_total)}], 'Location','best','FontSize',7);
grid on;
if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end
    fig1b_path = fullfile(fig_out_dir, 'algorithm1_psi_consensus.png');
    exportgraphics(gcf, fig1b_path, 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fig1b_path);
end

% Plot A3: Estimated n as xi evolves
% Note: xi_i(0)=0 for operators, so 1/xi_i(0)-1 = Inf at k=1.
% Replace non-finite values with NaN so MATLAB skips those points.
n_hat_hist = 1./xi_hist(1:n,:) - 1;

% Plot A_comp: Comparative results — four individual figures, one per metric
strat_labels = {'Coalition', 'RN-NE', 'Pareto-NE', 'Algo1'};
strat_cats   = categorical(strat_labels, strat_labels);  % preserve order
alpha_all = [alpha_coal, alpha_RN, alpha_sym, alpha_dist];
Pi_all    = [Pi_coal,    Pi_RN,    Pi_sym,    Pi_dist];
f_all     = [f_coal, f_RN, f_sym, f_dist];
g_all     = [g_coal, f_RN*(1-f_RN)*xr0^2, f_sym*(1-f_sym)*xr0^2, g_dist];

% A_comp_1: Bidding strategies alpha_i
figure('Position', [50, 50, 620, 460]);
bar(1:n, alpha_all);
xlabel('Operator'); ylabel('\alpha_i^*  (bid fraction)');
title('Bidding Strategies: All Four Methods');
xticklabels(op_labels);
yline(1, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
legend(strat_labels, 'Location', 'best', 'FontSize', 7);
grid on;
if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end
    exportgraphics(gcf, fullfile(fig_out_dir, 'algorithm1_comp_alpha_s1.png'), 'Resolution', 150);
    fprintf('\nFigure saved to: %s\n', fullfile(fig_out_dir, 'algorithm1_comp_alpha_s1.png'));
end

% A_comp_2: Individual profits Pi_i (S1)
figure('Position', [80, 80, 620, 460]);
bar(1:n, Pi_all);
xlabel('Operator'); ylabel('\Pi_i  ($/hr)');
title('Individual Profit: All Four Methods (S1)');
xticklabels(op_labels);
legend(strat_labels, 'Location', 'best', 'FontSize', 7);
grid on;
if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end
    exportgraphics(gcf, fullfile(fig_out_dir, 'algorithm1_comp_profit_s1.png'), 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fullfile(fig_out_dir, 'algorithm1_comp_profit_s1.png'));
end

% A_comp_2b: Per-unit profit Pi_i / x_i^o (S1)
Pipu_all_s1 = Pi_all ./ x0;   % $/hr per MW of installed capacity
figure('Position', [95, 95, 620, 460]);
bar(1:n, Pipu_all_s1);
xlabel('Operator'); ylabel('\Pi_i / x_i^o  ($/hr/MW)');
title('Per-Unit Profit: All Four Methods (S1, x^o=[120,100,80] MW)');
xticklabels(op_labels);
legend(strat_labels, 'Location', 'best', 'FontSize', 7);
grid on;
if save_figs
    style_current_figure();
    exportgraphics(gcf, fullfile(fig_out_dir, 'algorithm1_comp_profit_perunit_s1.png'), 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fullfile(fig_out_dir, 'algorithm1_comp_profit_perunit_s1.png'));
end


% =========================================================================
%  OVERBIDDING COMPARISON PLOTS — Scenario 2 params (x3^o = 30 MW)
%  Under S2: eta_3 = 0.12 < 1/(2n) = 1/6, so Pareto-NE and RN-NE both
%  overbid. Algorithm 1 detects saturation (if-else branch) and clamps
%  Op3 to alpha_3*=1 while re-optimising beta' for the remaining operators
%  via sigma/mu consensus (Corollary 5b hetero-beta), recovering f=0.5 exactly.
% =========================================================================

% --- Coalition under S2 ---
alpha_coal_s2 = 0.5 * ones(n, 1);
Pi_coal_s2    = profit(alpha_coal_s2, a2, a1, x0_s2, xL0, xr0_s2, eta_s2, 0,0,0);

% --- RN-NE under S2 (unclamped) ---
% nash_eq() clamps alpha, so compute the unclamped RN bid directly:
%   f_RN = n/(n+1) = 0.75,  alpha_i* = (1-f_RN)/eta_i = 0.25/eta_i
f_RN_s2_unc      = n / (n + 1);                          % = 0.75
alpha_RN_s2_unc  = (1 - f_RN_s2_unc) ./ eta_s2;         % [0.52, 0.63, 2.08]
Pi_RN_s2_unc     = profit(alpha_RN_s2_unc, a2, a1, x0_s2, xL0, xr0_s2, eta_s2, 0,0,0);

% --- Pareto-NE under S2 (unclamped symmetric Pareto, beta'=-1/6) ---
%     eta_3=0.12 => alpha_3* = 1/(6*0.12) = 1.39 > 1 => overbids
beta_sym_s2    = (1-n)/(4*n);                                    % = -1/6
f_sym_s2_unc   = (n + 4*n*beta_sym_s2)/(1 + n + 4*n*beta_sym_s2); % = 0.5 analytically
alpha_sym_s2_unc = (1 + 4*beta_sym_s2) * (1-f_sym_s2_unc) ./ eta_s2; % unclamped
Pi_sym_s2_unc  = profit(alpha_sym_s2_unc, a2, a1, x0_s2, xL0, xr0_s2, eta_s2, 0,0,0);

% --- Algo1 under S2 (Cor5b hetero-beta saturation branch) ---
%     alpha_opt_s2d = [0.396, 0.475, 1.000]: Op3 in Omega^c (alpha=1), f=0.5
Pi_dist_s2 = profit(alpha_opt_s2d, a2, a1, x0_s2, xL0, xr0_s2, eta_s2, 0,0,0);

% --- Unclamped OB under S2 (heterogeneous beta_cV, S2 params) ---
alpha_cV_unc_s2 = 0.5 * (1 + 4*beta_cV) ./ eta_s2;
Pi_cV_unc_s2    = profit(alpha_cV_unc_s2, a2, a1, x0_s2, xL0, xr0_s2, eta_s2, 0,0,0);

% --- No-overbid under S2 (clamp unclamped OB to [0,1]) ---
alpha_cV_clp_s2 = max(0, min(1, alpha_cV_unc_s2));
Pi_cV_clp_s2    = profit(alpha_cV_clp_s2, a2, a1, x0_s2, xL0, xr0_s2, eta_s2, 0,0,0);

% Assemble combined arrays — all evaluated at Scenario 2 params
ob_strat_labels = {'Coalition', 'RN-NE', 'Pareto-NE', 'Algo1 (Cor.5b)', 'Unclamped OB', 'No-overbid'};
n_ob = numel(ob_strat_labels);
ob_colors = lines(n_ob);
op_labels_ob = {sprintf('Op1 (%d MW)', x0_s2(1)), ...
                sprintf('Op2 (%d MW)', x0_s2(2)), ...
                sprintf('Op3 (%d MW)', x0_s2(3))};

alpha_combined = [alpha_coal_s2, alpha_RN_s2_unc, alpha_sym_s2_unc, alpha_opt_s2d, ...
                  alpha_cV_unc_s2, alpha_cV_clp_s2];               % n×6
Pi_combined    = [Pi_coal_s2, Pi_RN_s2_unc, Pi_sym_s2_unc, Pi_dist_s2, ...
                  Pi_cV_unc_s2, Pi_cV_clp_s2];                     % n×6
Pipu_combined  = Pi_combined ./ x0_s2;                             % per-unit ($/hr/MW)

% Plot A5: Bidding strategies — Scenario 2, all six strategies
% Pareto-NE and RN-NE bars for Op3 cross the alpha=1 boundary;
% Algo1 lands exactly at 1 (saturation branch activated).
figure('Position', [170, 170, 720, 480]);
b_al = bar(1:n, alpha_combined, 0.85);
for s = 1:n_ob, b_al(s).FaceColor = ob_colors(s,:); end
yline(1, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
% Cap y-axis at 1.6 for readability; annotate clipped RN-NE bar
ylim([0, 1.65]);
% Annotate Op3 bars that exceed the clipping limit
alpha_op3 = alpha_combined(3, :);
for s = 1:n_ob
    if alpha_op3(s) > 1.55
        text(3 + (s - (n_ob+1)/2)*0.14, 1.58, sprintf('%.2f', alpha_op3(s)), ...
            'HorizontalAlignment','center', 'FontSize', 6.5, 'Color', ob_colors(s,:), ...
            'FontWeight','bold');
    end
end
xlabel('Operator'); ylabel('\alpha_i^*  (bid fraction)');
title({'Bidding Strategies: Scenario 2 (x_3^o = 30 MW)'; ...
       'Algo.1 saturates Op3 \rightarrow \alpha_3^* = 1 (Corollary 5b)'});
xticks(1:n); xticklabels(op_labels_ob);
legend(ob_strat_labels, 'Location', 'northwest', 'FontSize', 7);
grid on;
if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end
    exportgraphics(gcf, fullfile(fig_out_dir, 'algorithm1_overbid_alpha_s2.png'), 'Resolution', 150);
    fprintf('\nFigure saved to: %s\n', fullfile(fig_out_dir, 'algorithm1_overbid_alpha_s2.png'));
end

% Plot A6: Per-operator absolute profit — Scenario 2, all six strategies
figure('Position', [200, 200, 720, 460]);
b_abs = bar(1:n, Pi_combined, 0.85);
for s = 1:n_ob, b_abs(s).FaceColor = ob_colors(s,:); end
xlabel('Operator'); ylabel('\Pi_i  ($/hr)');
title('Per-Operator Profit: Scenario 2 (x_3^o = 30 MW)');
xticks(1:n); xticklabels(op_labels_ob);
legend(ob_strat_labels, 'Location', 'northwest', 'FontSize', 7);
grid on;
if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end
    exportgraphics(gcf, fullfile(fig_out_dir, 'algorithm1_overbid_profit_s2.png'), 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fullfile(fig_out_dir, 'algorithm1_overbid_profit_s2.png'));
end

% Plot A7: Per-unit profit — Scenario 2, all six strategies
figure('Position', [230, 230, 720, 460]);
b_pu = bar(1:n, Pipu_combined, 0.85);
for s = 1:n_ob, b_pu(s).FaceColor = ob_colors(s,:); end
xlabel('Operator'); ylabel('\Pi_i / x_i^o  ($/hr/MW)');
title('Per-Unit Profit: Scenario 2 (x_3^o = 30 MW)');
xticks(1:n); xticklabels(op_labels_ob);
legend(ob_strat_labels, 'Location', 'northwest', 'FontSize', 7);
grid on;
if save_figs
    style_current_figure();
    exportgraphics(gcf, fullfile(fig_out_dir, 'algorithm1_overbid_profit_perunit_s2.png'), 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fullfile(fig_out_dir, 'algorithm1_overbid_profit_perunit_s2.png'));
end

% Visualization for Algorithm 2 (heterogeneous beta_i'(0))
if run_algo2
figure('Position', [80, 80, 1280, 820]);
sgtitle('Algorithm 2: Non-Equal \beta_i''(0) Coordination + One-Shot Bid', 'FontSize', 13);

subplot(2,3,1);
plot(1:beta_iters, Bsum_hist, 'LineWidth', 1.8); hold on;
yline(B_target, 'r--', 'LineWidth', 1.5);
xlabel('Outer iteration k'); ylabel('\Sigma_i \beta_i''(k)');
title('Aggregate Risk Parameter Convergence');
legend({'Estimated \Sigma\beta''', sprintf('Target=(1-n)/4=%.3f', B_target)}, 'Location','best');
grid on;

subplot(2,3,2);
plot(1:beta_iters, beta_hist', 'LineWidth', 1.5);
xlabel('Outer iteration k'); ylabel('\beta_i''(k)');
title('Individual \beta_i'' Trajectories');
legend(op_labels, 'Location','best'); grid on;

subplot(2,3,3);
plot(1:beta_iters, f2_hist, 'b', 'LineWidth', 1.8); hold on;
plot(1:beta_iters, g2_hist, 'm', 'LineWidth', 1.8);
yline(0.5, 'b--', 'LineWidth', 1.2);
yline(0.25*xr0^2, 'm--', 'LineWidth', 1.2);
xlabel('Outer iteration k'); ylabel('f(k), g(k)');
title('Efficiency Tracking');
legend({'f(k)', 'g(k)', 'f^*=0.5', 'g_{max}'}, 'Location','best');
grid on;

subplot(2,3,4);
plot(1:beta_iters, beta_err_hist, 'LineWidth', 1.8);
yline(0, 'k--', 'LineWidth', 1.3);
xlabel('Outer iteration k'); ylabel('e_\beta(k)');
title('Sum Error: e_\beta(k)=\Sigma\beta_i''(k)-\Sigma\beta_i''^*');
grid on;

subplot(2,3,5);
plot(1:beta_iters, sat2_hist, 'r', 'LineWidth', 1.8); hold on;
yline(0, 'k--', 'LineWidth', 1.2);
xlabel('Outer iteration k'); ylabel('count');
title('Saturation Count (\alpha_{unc}>1)');
ylim([-0.1, n+0.2]); grid on;

subplot(2,3,6);
plot(1:beta_iters, max(alpha2_unc_hist', [], 2), 'LineWidth', 1.8); hold on;
yline(1, 'k--', 'LineWidth', 1.2);
xlabel('Outer iteration k'); ylabel('max_i \alpha_i');
title('Peak Unclamped Bid');
grid on;

if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir')
        mkdir(fig_out_dir);
    end
    fig2_path = fullfile(fig_out_dir, 'algorithm2_results.png');
    exportgraphics(gcf, fig2_path, 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fig2_path);

    fig3_path = fullfile(fig_out_dir, 'algorithm2_bidding_comparison.png');
    fig3 = figure('Position', [120, 120, 1000, 460]);
    subplot(1,2,1);
    ba = bar(1:n, [alpha_RN, alpha_algo2]);
    ba(1).DisplayName = 'Risk-Neutral';
    ba(2).DisplayName = 'Algo 2 (unclamped)';
    xlabel('Operator'); ylabel('\alpha_i^*');
    title('Bidding Comparison');
    xticklabels(op_labels); legend('Location', 'best'); grid on;

    subplot(1,2,2);
    bp = bar(1:n, [Pi_RN, Pi_algo2]);
    bp(1).DisplayName = sprintf('Risk-Neutral (\\SigmaPi=%.0f)', sum(Pi_RN));
    bp(2).DisplayName = sprintf('Algo 2 unclamped (\\SigmaPi=%.0f)', sum(Pi_algo2));
    xlabel('Operator'); ylabel('\\Pi_i');
    title('Profit Comparison');
    xticklabels(op_labels); legend('Location', 'best'); grid on;
    style_current_figure();
    exportgraphics(fig3, fig3_path, 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fig3_path);

    fig4_path = fullfile(fig_out_dir, 'algorithm2_alpha_unclamped.png');
    fig4 = figure('Position', [120, 120, 1100, 700]);
    t = 1:beta_iters;
    for i = 1:n
        subplot(n,1,i);
        plot(t, alpha2_unc_hist(i,:), '-', 'LineWidth', 1.8); hold on;
        yline(1, 'k:', 'LineWidth', 1.2);
        ylabel(sprintf('\\alpha_%d', i));
        if i == 1
            title('Algorithm 2: Unclamped One-Shot Bids');
            legend({'unclamped', 'capacity limit'}, 'Location', 'best');
        end
        grid on;
    end
    xlabel('Outer iteration k');
    style_current_figure();
    exportgraphics(fig4, fig4_path, 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fig4_path);
end
end

% -------------------------------------------------------------------------
% Part 4 figures: Constrained RN vs Constrained ERS vs Cor5b hetero-beta (Scenario 2)
% -------------------------------------------------------------------------
op_labels_s2 = arrayfun(@(i) sprintf('Op%d',i), 1:n, 'UniformOutput', false);
scen4_labels = {'RN Unconstrained', 'Cor4 RN Clamped', 'Cor5b ERS Pareto'};
alpha_s4 = [alpha_RN_s2, alpha_Cor4, alpha_Rem1];
Pi_s4    = [Pi_RN_s2,    Pi_Cor4,   Pi_Rem1];
f_s4     = [f_RN_s2, f_Cor4, f_Rem1];
g_s4     = [g_RN_s2, g_Cor4, g_Rem1];

figure('Position', [50,50,680,460]);
bar(1:n, alpha_s4);
hold on; yline(1,'k--','LineWidth',1.5,'HandleVisibility','off');
xlabel('Operator'); ylabel('\alpha_i^*  (bid fraction)');
title('Corollary 4: Constrained RN vs ERS Pareto  (x_3^o=30 MW)');
xticklabels(op_labels_s2); legend(scen4_labels,'Location','best','FontSize',7);
grid on;
if save_figs
    style_current_figure();
    exportgraphics(gcf, fullfile(fig_out_dir, 'cor4_constrained_rn_alpha.png'), 'Resolution', 150);
    fprintf('Figure saved: cor4_constrained_rn_alpha.png\n');
end

% -------------------------------------------------------------------------
% Part 5 figures: (removed — Theorem 2 heterogeneous allocations not in paper)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Part 6 figures: Scenario 3 — six-strategy overbidding comparison
%   (analogous to Scenario 2 plots; Omega^c = {2,3})
% -------------------------------------------------------------------------
op_labels_s6 = arrayfun(@(i) sprintf('Op%d',i), 1:n, 'UniformOutput', false);

% --- Coalition under S3 ---
alpha_coal_s3 = 0.5 * ones(n, 1);
Pi_coal_s3    = profit(alpha_coal_s3, a2, a1, x0_s6, xL0, xr0_s6, eta_s6, 0,0,0);

% --- RN-NE under S3 (unclamped): f=n/(n+1)=0.75, alpha_i=0.25/eta_i ---
f_RN_s3_unc     = n / (n + 1);
alpha_RN_s3_unc = (1 - f_RN_s3_unc) ./ eta_s6;     % [0.357, 1.667, 1.667]
Pi_RN_s3_unc    = profit(alpha_RN_s3_unc, a2, a1, x0_s6, xL0, xr0_s6, eta_s6, 0,0,0);

% --- Pareto-NE under S3 (unclamped symmetric, already in alpha_unc_s6) ---
Pi_sym_s3_unc   = profit(alpha_unc_s6, a2, a1, x0_s6, xL0, xr0_s6, eta_s6, 0,0,0);

% --- Algo1 / Cor5b under S3 (alpha_opt_s6d = [0.286, 1, 1]) ---
Pi_dist_s3      = profit(alpha_opt_s6d, a2, a1, x0_s6, xL0, xr0_s6, eta_s6, 0,0,0);

% --- Unclamped OB under S3: heterogeneous beta, same as S2 deviator ---
beta_cV_s3      = [-0.20; -0.15; -0.05];
alpha_cV_unc_s3 = 0.5 * (1 + 4*beta_cV_s3) ./ eta_s6;
Pi_cV_unc_s3    = profit(alpha_cV_unc_s3, a2, a1, x0_s6, xL0, xr0_s6, eta_s6, 0,0,0);

% --- No-overbid under S3: naive per-operator clamp of Pareto-NE ---
alpha_cV_clp_s3 = max(0, min(1, alpha_unc_s6));
Pi_cV_clp_s3    = profit(alpha_cV_clp_s3, a2, a1, x0_s6, xL0, xr0_s6, eta_s6, 0,0,0);

ob_strat_labels_s3 = {'Coalition','RN-NE','Pareto-NE','Algo1 (Cor.5b)','Unclamped OB','No-overbid'};
n_ob_s3 = numel(ob_strat_labels_s3);
ob_colors_s3 = lines(n_ob_s3);
op_labels_s3 = {sprintf('Op1 (%d MW)', x0_s6(1)), ...
                sprintf('Op2 (%d MW)', x0_s6(2)), ...
                sprintf('Op3 (%d MW)', x0_s6(3))};

alpha_combined_s3 = [alpha_coal_s3, alpha_RN_s3_unc, alpha_unc_s6, alpha_opt_s6d, ...
                     alpha_cV_unc_s3, alpha_cV_clp_s3];
Pi_combined_s3    = [Pi_coal_s3, Pi_RN_s3_unc, Pi_sym_s3_unc, Pi_dist_s3, ...
                     Pi_cV_unc_s3, Pi_cV_clp_s3];
Pipu_combined_s3  = Pi_combined_s3 ./ x0_s6;

% Plot S3-A: Bidding strategies — all six strategies
figure('Position', [170, 170, 720, 480]);
b_al3 = bar(1:n, alpha_combined_s3, 0.85);
for s = 1:n_ob_s3, b_al3(s).FaceColor = ob_colors_s3(s,:); end
hold on; yline(1, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
ylim([0, 2.0]);
% Annotate bars that exceed axis clip
for s = 1:n_ob_s3
    for op = 1:n
        if alpha_combined_s3(op,s) > 1.85
            text(op + (s-(n_ob_s3+1)/2)*0.14, 1.88, ...
                sprintf('%.2f', alpha_combined_s3(op,s)), ...
                'HorizontalAlignment','center','FontSize',6.5, ...
                'Color',ob_colors_s3(s,:),'FontWeight','bold');
        end
    end
end
xlabel('Operator'); ylabel('\alpha_i^*  (bid fraction)');
title({'Bidding Strategies: Scenario 3  (x^o=[210,45,45] MW)', ...
       'Algo.1 saturates Ops 2,3 \rightarrow \alpha_{2,3}^* = 1  (Corollary 5b)'});
xticks(1:n); xticklabels(op_labels_s3);
legend(ob_strat_labels_s3, 'Location', 'northwest', 'FontSize', 7);
grid on;
if save_figs
    style_current_figure();
    exportgraphics(gcf, fullfile(fig_out_dir, 'algorithm1_overbid_alpha_s3.png'), 'Resolution', 150);
    fprintf('Figure saved: algorithm1_overbid_alpha_s3.png\n');
end

% Plot S3-B: Per-operator profit — all six strategies
figure('Position', [200, 200, 720, 460]);
b_abs3 = bar(1:n, Pi_combined_s3, 0.85);
for s = 1:n_ob_s3, b_abs3(s).FaceColor = ob_colors_s3(s,:); end
xlabel('Operator'); ylabel('\Pi_i  ($/hr)');
title('Per-Operator Profit: Scenario 3  (x^o=[210,45,45] MW)');
xticks(1:n); xticklabels(op_labels_s3);
legend(ob_strat_labels_s3, 'Location', 'northwest', 'FontSize', 7);
grid on;
if save_figs
    style_current_figure();
    exportgraphics(gcf, fullfile(fig_out_dir, 'algorithm1_overbid_profit_s3.png'), 'Resolution', 150);
    fprintf('Figure saved: algorithm1_overbid_profit_s3.png\n');
end

% Plot S3-C: Per-unit profit — all six strategies
figure('Position', [230, 230, 720, 460]);
b_pu3 = bar(1:n, Pipu_combined_s3, 0.85);
for s = 1:n_ob_s3, b_pu3(s).FaceColor = ob_colors_s3(s,:); end
xlabel('Operator'); ylabel('\Pi_i / x_i^o  ($/hr/MW)');
title('Per-Unit Profit: Scenario 3  (x^o=[210,45,45] MW)');
xticks(1:n); xticklabels(op_labels_s3);
legend(ob_strat_labels_s3, 'Location', 'northwest', 'FontSize', 7);
grid on;
if save_figs
    style_current_figure();
    exportgraphics(gcf, fullfile(fig_out_dir, 'algorithm1_overbid_profit_perunit_s3.png'), 'Resolution', 150);
    fprintf('Figure saved: algorithm1_overbid_profit_perunit_s3.png\n');
end

% -------------------------------------------------------------------------
% Part 7 figures: Phase 4 xi2/psi2 consensus trajectories (Cor. 5b branch)
% -------------------------------------------------------------------------
% Note: Phase 4 uses the same xi/psi ratio trick as Phase 1, seeded only by
% operators NOT in Omega_ERS.  Figures are omitted here; the xi/psi consensus
% plots for Phase 1 (above) already illustrate the protocol mechanics.


if close_figs
    close all;
    fprintf('\nAll figures closed (close_figs=true).\n');
end

%% ========================================================================
%                         HELPER FUNCTIONS
% =========================================================================

function [alpha_star, f, Pi, J] = nash_eq(beta_prime, n, a2, a1, x0, xL0, xr0, eta, sigma_r2, sigma_L2, sigma_rL)
    % Nash equilibrium given beta' using closed-form from Theorem 2 proof:
    %   f* = (n + 4*sum(beta')) / (1 + n + 4*sum(beta'))
    %   alpha_i* = (1 + 4*beta_i') * (1 - f*) / eta_i
    %
    % This matches eq. alpha_equation after substituting the Pareto condition.
    B_sum = sum(beta_prime);
    f = (n + 4*B_sum) / (1 + n + 4*B_sum);
    alpha_star = (1 + 4*beta_prime) .* (1-f) ./ eta;
    alpha_star = max(0, min(1, alpha_star));
    f = sum(alpha_star .* eta);   % recompute f from (possibly clamped) alpha
    Pi = profit(alpha_star, a2, a1, x0, xL0, xr0, eta, sigma_r2, sigma_L2, sigma_rL);
    price_var = 4*a2^2 * (xr0^2*(1-f)^2 + sigma_L2 + sigma_r2 - 2*sigma_rL);
    J = Pi - (beta_prime/a2) .* price_var;
end

function Pi = profit(alpha, a2, a1, x0, xL0, xr0, eta, sigma_r2, sigma_L2, sigma_rL)
    % Individual profits (eq. Pi_exp):
    %   Pi_i = 2*a2*alpha_i*x_i^o*x_r^o*(1-f) + a1*x_i^o + 2*a2*x_i^o*(xL0-xr0)
    f  = sum(alpha .* eta);
    Pi = 2*a2 * alpha .* x0 * xr0 .* (1-f) ...
       + a1 * x0 ...
       + 2*a2 * x0 .* (xL0 - xr0);
    if sigma_rL ~= 0 || sigma_r2 > 0
        Pi = Pi + 2*a2*(sigma_rL - sigma_r2/length(x0)) * ones(length(x0),1);
    end
end

function style_current_figure()
    set(gcf, 'Color', 'w', 'InvertHardcopy', 'off');
    ax = findall(gcf, 'Type', 'axes');
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    for i = 1:numel(ax)
        t = get(ax(i), 'Title');
        set(t, 'Color', 'k');
        xl = get(ax(i), 'XLabel');
        yl = get(ax(i), 'YLabel');
        set(xl, 'Color', 'k');
        set(yl, 'Color', 'k');
    end
    lgd = findall(gcf, 'Type', 'legend');
    set(lgd, 'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'k');
    tx = findall(gcf, 'Type', 'text');
    set(tx, 'Color', 'k');
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

function [z_final, z_hist, n_iters] = run_cons(z_init, L, eps_c, max_iter, tol)
    % Run a consensus protocol z(k+1) = z(k) - eps*L*z(k) to convergence.
    z = z_init;  z_hist = zeros(length(z_init), max_iter);  n_iters = max_iter;
    for k = 1:max_iter
        z_hist(:,k) = z;
        z_new = z - eps_c * L * z;
        if k>1 && norm(z_new-z) < tol,  n_iters=k;  z=z_new;  break;  end
        z = z_new;
    end
    z_final = z;
end
