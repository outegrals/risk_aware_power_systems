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

save_figs = true;
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
%              PART 3: ALGORITHM 1 OVERBIDDING SENSITIVITY
% =========================================================================

fprintf('\n----------------------------------------------------------\n');
fprintf('  PART 3: ALGORITHM 1 OVERBIDDING SENSITIVITY\n');
fprintf('----------------------------------------------------------\n\n');

Nb = length(beta_sweep);
alpha_unc_hist_b = zeros(n, Nb);
alpha_clp_hist_b = zeros(n, Nb);
f_unc_b = zeros(1, Nb);
f_clp_b = zeros(1, Nb);
g_unc_b = zeros(1, Nb);
g_clp_b = zeros(1, Nb);
Pi_unc_b = zeros(1, Nb);
Pi_clp_b = zeros(1, Nb);
Pi_unc_i_b = zeros(n, Nb);
Pi_clp_i_b = zeros(n, Nb);
sat_b = zeros(1, Nb);

for k = 1:Nb
    bp = beta_sweep(k);
    alpha_unc = 0.5 * (1 + 4*bp) ./ eta_dist;
    alpha_clp = max(0, min(1, alpha_unc));

    alpha_unc_hist_b(:,k) = alpha_unc;
    alpha_clp_hist_b(:,k) = alpha_clp;
    sat_b(k) = sum(alpha_unc > 1);

    f_unc_b(k) = sum(alpha_unc .* eta_dist);
    f_clp_b(k) = sum(alpha_clp .* eta_dist);
    g_unc_b(k) = f_unc_b(k) * (1 - f_unc_b(k)) * xr0^2;
    g_clp_b(k) = f_clp_b(k) * (1 - f_clp_b(k)) * xr0^2;

    Pi_unc_i = profit(alpha_unc, a2, a1, x0, xL0, xr0, eta_dist, sigma_r2, sigma_L2, sigma_rL);
    Pi_clp_i = profit(alpha_clp, a2, a1, x0, xL0, xr0, eta_dist, sigma_r2, sigma_L2, sigma_rL);
    Pi_unc_i_b(:,k) = Pi_unc_i;
    Pi_clp_i_b(:,k) = Pi_clp_i;
    Pi_unc_b(k) = sum(Pi_unc_i);
    Pi_clp_b(k) = sum(Pi_clp_i);
end

[Pi_unc_max, idx_unc_max] = max(Pi_unc_b);
[Pi_clp_max, idx_clp_max] = max(Pi_clp_b);
fprintf('Unclamped max sum(Pi)=%.2f at beta''=%.4f\n', Pi_unc_max, beta_sweep(idx_unc_max));
fprintf('No-overbidding max sum(Pi)=%.2f at beta''=%.4f\n', Pi_clp_max, beta_sweep(idx_clp_max));
fprintf('Pareto-central reference sum(Pi)=%.2f\n', sum(Pi_sym));

% Per-operator overbidding thresholds (eq:overbid_bound): alpha_i > 1 iff beta_i' > (2*eta_i-1)/4
overbid_thresh = (2*eta_dist - 1) / 4;
fprintf('\nPer-operator saturation thresholds (alpha_i > 1 iff beta'' > threshold):\n');
for i = 1:n
    fprintf('  Op %d: eta=%.4f => threshold beta'' = %.4f\n', i, eta_dist(i), overbid_thresh(i));
end

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
fprintf('\n  Overbidding sweep (Algorithm 1):\n');
fprintf('    max sum(Pi), unclamped          = %.2f at beta''=%.4f\n', Pi_unc_max, beta_sweep(idx_unc_max));
fprintf('    max sum(Pi), no-overbidding     = %.2f at beta''=%.4f\n', Pi_clp_max, beta_sweep(idx_clp_max));

%% ========================================================================
%                         VISUALIZATION
% =========================================================================

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
xlabel('Iteration k'); ylabel('\psi_i(k)');
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
n_hat_plot = n_hat_hist(:,1:K);
n_hat_plot(~isfinite(n_hat_plot)) = NaN;
figure('Position', [110, 110, 620, 460]);
plot(1:K, n_hat_plot', 'LineWidth', 1.5);
yline(n, 'k--', 'LineWidth', 1.5);
xlabel('Iteration k'); ylabel('\hat{n}_i(k)');
title('Distributed Estimate of n');
legend([op_labels, {sprintf('True n=%d',n)}], 'Location','best','FontSize',7);
ylim([0, n + 3]); grid on;   % fixed ylim: transient overshoots n, but converges to n
if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end
    fig1c_path = fullfile(fig_out_dir, 'algorithm1_nhat_estimate.png');
    exportgraphics(gcf, fig1c_path, 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fig1c_path);
end

% Plot A4: Estimated x_r^o as ratio psi/xi evolves
% Note: xi_i(0)=0, so psi_i(0)/xi_i(0) = Inf at k=1; early iterates are very large.
% Replace non-finite values with NaN and fix ylim so convergence is visible.
xr0_hist_ops = psi_hist(1:n,:) ./ xi_hist(1:n,:);
xr0_hist_plot = xr0_hist_ops(:,1:K);
xr0_hist_plot(~isfinite(xr0_hist_plot)) = NaN;
figure('Position', [140, 140, 620, 460]);
plot(1:K, xr0_hist_plot', 'LineWidth', 1.5);
yline(xr0, 'k--', 'LineWidth', 1.5);
xlabel('Iteration k'); ylabel('\hat{x}_r^o estimate');
title('Distributed Estimate of x_r^o');
legend([op_labels, {sprintf('True=%.0f',xr0)}], 'Location','best','FontSize',7);
ylim([0, xr0 * 2]); grid on;  % fixed ylim: clips large transient, shows convergence
if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end
    fig1d_path = fullfile(fig_out_dir, 'algorithm1_xr0_estimate.png');
    exportgraphics(gcf, fig1d_path, 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fig1d_path);
end

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
xlabel('Operator'); ylabel('\alpha_i^*');
title('Bidding Strategies: All Four Methods');
xticklabels(op_labels);
yline(1, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
legend(strat_labels, 'Location', 'best', 'FontSize', 7);
grid on;
if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end
    exportgraphics(gcf, fullfile(fig_out_dir, 'algorithm1_comp_alpha.png'), 'Resolution', 150);
    fprintf('\nFigure saved to: %s\n', fullfile(fig_out_dir, 'algorithm1_comp_alpha.png'));
end

% A_comp_2: Individual profits Pi_i
figure('Position', [80, 80, 620, 460]);
bar(1:n, Pi_all);
xlabel('Operator'); ylabel('\Pi_i');
title('Individual Profit: All Four Methods');
xticklabels(op_labels);
legend(strat_labels, 'Location', 'best', 'FontSize', 7);
grid on;
if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end
    exportgraphics(gcf, fullfile(fig_out_dir, 'algorithm1_comp_profit.png'), 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fullfile(fig_out_dir, 'algorithm1_comp_profit.png'));
end

% A_comp_3: Aggregate bid fraction f
figure('Position', [110, 110, 620, 460]);
bar(strat_cats, f_all);
yline(0.5, 'r--', 'LineWidth', 1.5);
xlabel('Strategy'); ylabel('f(\alpha)');
title('Aggregate Bid Fraction (target f = 0.5)');
ylim([0, 1]); grid on;
if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end
    exportgraphics(gcf, fullfile(fig_out_dir, 'algorithm1_comp_f.png'), 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fullfile(fig_out_dir, 'algorithm1_comp_f.png'));
end

% A_comp_4: Market efficiency g
figure('Position', [140, 140, 620, 460]);
bar(strat_cats, g_all);
yline(0.25*xr0^2, 'r--', 'LineWidth', 1.5);
xlabel('Strategy'); ylabel('g(\alpha)  [MW^2]');
title(sprintf('Market Efficiency g  (g_{max} = %.0f MW^2)', 0.25*xr0^2));
grid on;
if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end
    exportgraphics(gcf, fullfile(fig_out_dir, 'algorithm1_comp_g.png'), 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fullfile(fig_out_dir, 'algorithm1_comp_g.png'));
end

% Plot A5: overbidding bids
figure('Position', [170, 170, 700, 460]);
plot(beta_sweep, alpha_unc_hist_b', 'LineWidth', 1.4); hold on;
yline(1, 'k--', 'LineWidth', 1.2);
xlabel('\beta'' (common)'); ylabel('\alpha_i (unclamped)');
title('Unclamped Bids vs \beta''');
legend([op_labels, {'Capacity limit'}], 'Location', 'best');
grid on;
if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end
    fig1e_path = fullfile(fig_out_dir, 'algorithm1_overbid_alpha.png');
    exportgraphics(gcf, fig1e_path, 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fig1e_path);
end

% Plot A6: overbidding count
figure('Position', [200, 200, 700, 460]);
plot(beta_sweep, sat_b, 'r', 'LineWidth', 1.8); hold on;
yline(0, 'k--', 'LineWidth', 1.1);
xlabel('\beta'' (common)'); ylabel('count');
title('Number of Overbidding Operators');
ylim([-0.1, n+0.2]); grid on;

if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir')
        mkdir(fig_out_dir);
    end
    fig1f_path = fullfile(fig_out_dir, 'algorithm1_overbid_count.png');
    exportgraphics(gcf, fig1f_path, 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fig1f_path);
end

% Plot A7: overbidding profit impact
figure('Position', [230, 230, 700, 460]);
plot(beta_sweep, Pi_unc_b, 'b', 'LineWidth', 1.8); hold on;
plot(beta_sweep, Pi_clp_b, 'm', 'LineWidth', 1.8);
yline(sum(Pi_RN), 'k:', 'LineWidth', 1.2);
yline(sum(Pi_sym), 'g--', 'LineWidth', 1.2);
xlabel('\beta'' (common)'); ylabel('\Sigma \Pi_i');
title('Profit Impact: Unclamped vs No-Overbidding');
legend({'Unclamped', 'No-overbidding', 'RN ref', 'Pareto ref'}, 'Location', 'best');
grid on;

if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir')
        mkdir(fig_out_dir);
    end
    fig1g_path = fullfile(fig_out_dir, 'algorithm1_overbid_profit.png');
    exportgraphics(gcf, fig1g_path, 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fig1g_path);
end

% Plot A8: per-operator overbidding profit impact
figure('Position', [250, 250, 760, 860]);
tiledlayout(n, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:n
    nexttile;
    plot(beta_sweep, Pi_unc_i_b(i,:), 'b', 'LineWidth', 1.6); hold on;
    plot(beta_sweep, Pi_clp_i_b(i,:), 'm', 'LineWidth', 1.6);
    yline(Pi_RN(i), 'k:', 'LineWidth', 1.1);
    yline(Pi_sym(i), 'g--', 'LineWidth', 1.1);
    xline(overbid_thresh(i), 'r:', 'LineWidth', 1.3);
    ylabel(sprintf('\\Pi_%d', i));
    title(sprintf('Operator %d Profit vs \\beta'' (sat.thresh=%.3f)', i, overbid_thresh(i)));
    if i == 1
        legend({'Unclamped', 'No-overbidding', 'RN ref', 'Pareto ref', 'Sat. threshold'}, 'Location', 'best');
    end
    grid on;
end
xlabel('\beta'' (common)');
if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir')
        mkdir(fig_out_dir);
    end
    fig1h_path = fullfile(fig_out_dir, 'algorithm1_overbid_profit_individual.png');
    exportgraphics(gcf, fig1h_path, 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fig1h_path);
end

% Plot A9: overbidding aggregate f
figure('Position', [280, 280, 700, 460]);
plot(beta_sweep, f_unc_b, 'b', 'LineWidth', 1.6); hold on;
plot(beta_sweep, f_clp_b, 'm', 'LineWidth', 1.6);
yline(0.5, 'k--', 'LineWidth', 1.2);
xlabel('\beta'' (common)'); ylabel('f(\alpha)');
title('Aggregate Bid Fraction');
legend({'Unclamped', 'No-overbidding', 'Pareto target'}, 'Location', 'best');
grid on;

if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir')
        mkdir(fig_out_dir);
    end
    fig1i_path = fullfile(fig_out_dir, 'algorithm1_overbid_f.png');
    exportgraphics(gcf, fig1i_path, 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fig1i_path);
end

% Plot A10: Comprehensive per-operator overbidding effects (algorithm1_overbidding.png)
% Row 1: per-operator alpha_i (unclamped vs clamped) vs beta' with saturation threshold
% Row 2: per-operator Pi_i (unclamped vs clamped) vs beta' with saturation threshold
figure('Position', [50, 50, 1200, 820]);
tiledlayout(2, n, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:n
    nexttile;
    plot(beta_sweep, alpha_unc_hist_b(i,:), 'b', 'LineWidth', 1.6); hold on;
    plot(beta_sweep, alpha_clp_hist_b(i,:), 'm--', 'LineWidth', 1.6);
    xline(overbid_thresh(i), 'r:', 'LineWidth', 1.5);
    yline(1, 'k--', 'LineWidth', 1.2);
    xlabel('\beta'' (common)'); ylabel(sprintf('\\alpha_%d', i));
    title(sprintf('Op %d bid (thresh=%.3f)', i, overbid_thresh(i)));
    if i == 1
        legend({'Unclamped', 'No-overbid', 'Sat. threshold', 'Cap. limit'}, 'FontSize', 7, 'Location', 'best');
    end
    grid on;
end
for i = 1:n
    nexttile;
    plot(beta_sweep, Pi_unc_i_b(i,:), 'b', 'LineWidth', 1.6); hold on;
    plot(beta_sweep, Pi_clp_i_b(i,:), 'm--', 'LineWidth', 1.6);
    xline(overbid_thresh(i), 'r:', 'LineWidth', 1.5);
    yline(Pi_RN(i), 'k:', 'LineWidth', 1.1);
    yline(Pi_sym(i), 'g--', 'LineWidth', 1.1);
    xlabel('\beta'' (common)'); ylabel(sprintf('\\Pi_%d', i));
    title(sprintf('Op %d profit', i));
    if i == 1
        legend({'Unclamped', 'No-overbid', 'Sat. threshold', 'RN ref', 'Pareto ref'}, 'FontSize', 7, 'Location', 'best');
    end
    grid on;
end
if save_figs
    style_current_figure();
    if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end
    exportgraphics(gcf, fullfile(fig_out_dir, 'algorithm1_overbidding.png'), 'Resolution', 150);
    fprintf('Figure saved to: %s\n', fullfile(fig_out_dir, 'algorithm1_overbidding.png'));
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
