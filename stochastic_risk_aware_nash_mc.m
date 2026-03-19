%% Stochastic Risk-Aware Nash Study (Monte Carlo Reference)
% This file is intentionally separate from existing baseline scripts.
% It keeps the same analytical bidding formulas, but evaluates outcomes
% with sampled real-time realizations x_i and x_L so sigma terms matter.

clc; clear; close all;
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesColor', 'w');
set(groot, 'defaultAxesXColor', 'k');
set(groot, 'defaultAxesYColor', 'k');
set(groot, 'defaultTextColor', 'k');
set(groot, 'defaultLegendTextColor', 'k');
set(groot, 'defaultLegendColor', 'w');
set(groot, 'defaultLegendEdgeColor', 'k');

%% ========================================================================
%                            USER PARAMETERS
% =========================================================================

% Market parameters
n  = 3;
a2 = 0.1;
a1 = 1;

% Day-ahead forecasts
x0  = [120; 100; 80];   % renewable means
xL0 = 400;              % load mean
xr0 = sum(x0);
eta = x0 / xr0;

% Transformed risk parameters beta' = a2 * beta
beta_prime_RN  = zeros(n,1);
beta_prime_opt = ((1 - n) / (4*n)) * ones(n,1);

% Monte Carlo configuration
Nmc      = 2e5;       % number of scenarios
rng_seed = 11;        % reproducibility

% Uncertainty model for [x_1, ..., x_n, x_L]
% Choose std devs and correlation here (adjust as needed).
sigma_x = [18; 14; 10];   % std dev of renewable operators
sigma_L = 25;             % std dev of load
rho_rr  = 0.20;           % pairwise renewable-renewable correlation
rho_rL  = 0.35;           % renewable-load correlation

% Overbidding sweep (common beta')
beta_sweep = linspace(-0.24, 0.35, 180);

% Algorithm 1 (distributed one-shot consensus) settings
run_algo1      = true;
epsilon_c      = 0.1;
max_cons_iter  = 500;
cons_tol       = 1e-10;

save_figs = true;
fig_out_dir = 'figures';
if save_figs && ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

%% ========================================================================
%                         BUILD COVARIANCE MODEL
% =========================================================================

mu = [x0; xL0];
Sigma = zeros(n+1);

for i = 1:n
    Sigma(i,i) = sigma_x(i)^2;
end
Sigma(n+1,n+1) = sigma_L^2;

for i = 1:n
    for j = 1:n
        if i ~= j
            Sigma(i,j) = rho_rr * sigma_x(i) * sigma_x(j);
        end
    end
    Sigma(i,n+1) = rho_rL * sigma_x(i) * sigma_L;
    Sigma(n+1,i) = Sigma(i,n+1);
end

% PSD guard for numerical robustness.
Sigma = nearest_psd(Sigma);

%% ========================================================================
%                    ANALYTICAL ALPHAS (UNCHANGED THEORY)
% =========================================================================

[alpha_RN, f_RN]   = nash_alpha(beta_prime_RN, n, eta);
[alpha_opt, f_opt] = nash_alpha(beta_prime_opt, n, eta);

fprintf('============================================================\n');
fprintf('   STOCHASTIC MONTE CARLO CHECK (separate reference file)\n');
fprintf('============================================================\n');
fprintf('Nmc = %d, seed = %d\n', Nmc, rng_seed);
fprintf('RN alpha      = [%s], f=%.6f\n', num2str(alpha_RN',  '%.6f '), f_RN);
fprintf('Pareto alpha  = [%s], f=%.6f\n', num2str(alpha_opt', '%.6f '), f_opt);

if run_algo1
    algo1 = run_algorithm1(x0, n, epsilon_c, max_cons_iter, cons_tol);
    fprintf('\nAlgorithm 1 (distributed one-shot) estimates:\n');
    fprintf('  n_hat      = [%s]\n', num2str(algo1.n_hat', '%.6f '));
    fprintf('  eta_hat    = [%s]\n', num2str(algo1.eta_hat', '%.6f '));
    fprintf('  beta''_hat  = [%s]\n', num2str(algo1.beta_prime_hat', '%.6f '));
    fprintf('  alpha_unc  = [%s]\n', num2str(algo1.alpha_unc', '%.6f '));
    fprintf('  alpha_clp  = [%s]\n', num2str(algo1.alpha_clp', '%.6f '));
end

%% ========================================================================
%                 MONTE CARLO: RANDOM REAL-TIME REALIZATIONS
% =========================================================================

rng(rng_seed);
X = mvnrnd_chol(mu, Sigma, Nmc);   % Nmc x (n+1)
Xr = X(:,1:n);
XL = X(:,n+1);

% Prevent nonphysical negative generation/load in samples.
Xr = max(0, Xr);
XL = max(0, XL);

% Evaluate RN and Pareto profiles under sampled uncertainty.
stats_RN  = evaluate_profile(alpha_RN,  beta_prime_RN,  x0, xL0, a2, a1, Xr, XL);
stats_opt = evaluate_profile(alpha_opt, beta_prime_opt, x0, xL0, a2, a1, Xr, XL);
if run_algo1
    stats_algo1_unc = evaluate_profile(algo1.alpha_unc, algo1.beta_prime_hat, x0, xL0, a2, a1, Xr, XL);
    stats_algo1_clp = evaluate_profile(algo1.alpha_clp, algo1.beta_prime_hat, x0, xL0, a2, a1, Xr, XL);
end

fprintf('\nMonte Carlo expected outcomes:\n');
print_stats('Risk-neutral', stats_RN);
print_stats('Pareto-risk',  stats_opt);
if run_algo1
    print_stats('Algorithm1-unclamped', stats_algo1_unc);
    print_stats('Algorithm1-clamped',   stats_algo1_clp);
end

%% ========================================================================
%              OVERBIDDING SENSITIVITY WITH STOCHASTIC PAYOFFS
% =========================================================================

Nb = numel(beta_sweep);
sumPi_unc = zeros(Nb,1);
sumPi_clp = zeros(Nb,1);
Jsum_unc  = zeros(Nb,1);
Jsum_clp  = zeros(Nb,1);
sat_count = zeros(Nb,1);

for k = 1:Nb
    bp = beta_sweep(k);
    bpv = bp * ones(n,1);

    alpha_unc = 0.5 * (1 + 4*bpv) ./ eta;
    alpha_clp = max(0, min(1, alpha_unc));

    sat_count(k) = sum(alpha_unc > 1);

    st_unc = evaluate_profile(alpha_unc, bpv, x0, xL0, a2, a1, Xr, XL);
    st_clp = evaluate_profile(alpha_clp, bpv, x0, xL0, a2, a1, Xr, XL);

    sumPi_unc(k) = sum(st_unc.Pi_mean_i);
    sumPi_clp(k) = sum(st_clp.Pi_mean_i);
    Jsum_unc(k)  = sum(st_unc.J_i);
    Jsum_clp(k)  = sum(st_clp.J_i);
end

%% ========================================================================
%                                PLOTS
% =========================================================================

if run_algo1
    fig0 = figure('Color','w');
    tiledlayout(2,1);
    nexttile;
    plot(1:algo1.xi_iters, algo1.xi_hist(1:n,1:algo1.xi_iters)', 'LineWidth', 1.4); hold on;
    yline(1/(n+1), 'k--', 'LineWidth', 1.2);
    xlabel('Iteration k'); ylabel('\xi_i(k)');
    title('Algorithm 1 Consensus: \xi_i');
    grid on;
    nexttile;
    plot(1:algo1.psi_iters, algo1.psi_hist(1:n,1:algo1.psi_iters)', 'LineWidth', 1.4); hold on;
    yline(sum(x0)/(n+1), 'k--', 'LineWidth', 1.2);
    xlabel('Iteration k'); ylabel('\psi_i(k)');
    title('Algorithm 1 Consensus: \psi_i');
    grid on;
    style_current_figure();
    if save_figs
        exportgraphics(fig0, fullfile(fig_out_dir, 'mc_algorithm1_consensus.png'), 'Resolution', 150);
    end
end

fig1 = figure('Color','w');
plot(beta_sweep, sumPi_unc, 'LineWidth', 1.8); hold on;
plot(beta_sweep, sumPi_clp, 'LineWidth', 1.8);
xlabel('\beta'' (common)');
ylabel('E[\Sigma_i \Pi_i]');
title('Monte Carlo Expected Aggregate Profit');
legend({'unclamped', 'clamped [0,1]'}, 'Location', 'best');
grid on;
style_current_figure();
if save_figs
    exportgraphics(fig1, fullfile(fig_out_dir, 'mc_overbid_profit.png'), 'Resolution', 150);
end

fig2 = figure('Color','w');
plot(beta_sweep, Jsum_unc, 'LineWidth', 1.8); hold on;
plot(beta_sweep, Jsum_clp, 'LineWidth', 1.8);
xlabel('\beta'' (common)');
ylabel('\Sigma_i J_i');
title('Monte Carlo Risk-Aware Objective');
legend({'unclamped', 'clamped [0,1]'}, 'Location', 'best');
grid on;
style_current_figure();
if save_figs
    exportgraphics(fig2, fullfile(fig_out_dir, 'mc_overbid_Jsum.png'), 'Resolution', 150);
end

fig3 = figure('Color','w');
plot(beta_sweep, sat_count, 'r', 'LineWidth', 1.8);
xlabel('\beta'' (common)');
ylabel('count');
title('Overbidding Count (Unclamped \alpha_i > 1)');
grid on;
style_current_figure();
if save_figs
    exportgraphics(fig3, fullfile(fig_out_dir, 'mc_overbid_count.png'), 'Resolution', 150);
end

if save_figs
    fprintf('\nSaved MC figures to %s/\n', fig_out_dir);
end

%% ========================================================================
%                            HELPER FUNCTIONS
% =========================================================================

function [alpha_star, f] = nash_alpha(beta_prime, n, eta)
    B_sum = sum(beta_prime);
    f = (n + 4*B_sum) / (1 + n + 4*B_sum);
    alpha_star = (1 + 4*beta_prime) .* (1 - f) ./ eta;
end

function out = run_algorithm1(x0, n, epsilon_c, max_cons_iter, cons_tol)
    N_total = n + 1;
    Adj = ones(N_total) - eye(N_total);
    L = diag(sum(Adj,2)) - Adj;

    xi = zeros(N_total,1);   xi(n+1) = 1;
    psi = zeros(N_total,1);  psi(1:n) = x0;

    xi_hist = zeros(N_total, max_cons_iter);
    psi_hist = zeros(N_total, max_cons_iter);
    xi_iters = max_cons_iter;
    psi_iters = max_cons_iter;

    for k = 1:max_cons_iter
        xi_hist(:,k) = xi;
        psi_hist(:,k) = psi;
        xi_new = xi - epsilon_c * L * xi;
        psi_new = psi - epsilon_c * L * psi;

        if k > 1
            if norm(xi_new - xi) < cons_tol && xi_iters == max_cons_iter
                xi_iters = k;
            end
            if norm(psi_new - psi) < cons_tol && psi_iters == max_cons_iter
                psi_iters = k;
            end
        end
        xi = xi_new;
        psi = psi_new;
    end

    xi_f = xi(1:n);
    psi_f = psi(1:n);

    n_hat = 1 ./ xi_f - 1;
    xr0_hat = psi_f ./ xi_f;
    eta_hat = x0 ./ xr0_hat;
    beta_prime_hat = (1 - n_hat) ./ (4 * n_hat);

    alpha_unc = 0.5 * (1 + 4 * beta_prime_hat) ./ eta_hat;
    alpha_clp = max(0, min(1, alpha_unc));

    out.n_hat = n_hat;
    out.eta_hat = eta_hat;
    out.beta_prime_hat = beta_prime_hat;
    out.alpha_unc = alpha_unc;
    out.alpha_clp = alpha_clp;
    out.xi_hist = xi_hist;
    out.psi_hist = psi_hist;
    out.xi_iters = xi_iters;
    out.psi_iters = psi_iters;
end

function st = evaluate_profile(alpha, beta_prime, x0, xL0, a2, a1, Xr, XL)
    n = numel(alpha);
    xr0 = sum(x0);

    lambda0 = 2*a2*(xL0 - sum(alpha .* x0)) + a1;
    lambda  = 2*a2*(XL - sum(Xr,2)) + a1;

    Pi_samp = zeros(size(Xr));
    for i = 1:n
        Pi_samp(:,i) = lambda0 * alpha(i) * x0(i) + lambda .* (Xr(:,i) - alpha(i)*x0(i));
    end

    dlam2 = (lambda - lambda0).^2;
    Pi_mean_i = mean(Pi_samp, 1)';
    J_i = Pi_mean_i - (beta_prime / a2) * mean(dlam2);

    f = sum(alpha .* (x0 / xr0));
    g = f * (1 - f) * xr0^2;

    st.alpha = alpha;
    st.f = f;
    st.g = g;
    st.lambda0 = lambda0;
    st.lambda_mean = mean(lambda);
    st.lambda_var = var(lambda, 1);
    st.Pi_mean_i = Pi_mean_i;
    st.Pi_var_i  = var(Pi_samp, 1, 1)';
    st.J_i = J_i;
end

function X = mvnrnd_chol(mu, Sigma, N)
    d = numel(mu);
    L = chol(Sigma, 'lower');
    Z = randn(N, d);
    X = Z * L' + repmat(mu(:)', N, 1);
end

function A = nearest_psd(A)
    A = (A + A') / 2;
    [V,D] = eig(A);
    d = diag(D);
    d(d < 1e-12) = 1e-12;
    A = V * diag(d) * V';
    A = (A + A') / 2;
end

function print_stats(label, st)
    fprintf('  %s:\n', label);
    fprintf('    E[lambda] = %.4f, Var(lambda) = %.4f\n', st.lambda_mean, st.lambda_var);
    fprintf('    f(alpha) = %.6f, g(alpha) = %.2f\n', st.f, st.g);
    fprintf('    E[sum Pi] = %.2f, sum J = %.2f\n', sum(st.Pi_mean_i), sum(st.J_i));
end

function style_current_figure()
    set(gcf, 'Color', 'w', 'InvertHardcopy', 'off');
    ax = findall(gcf, 'Type', 'axes');
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    for i = 1:numel(ax)
        t = get(ax(i), 'Title'); set(t, 'Color', 'k');
        xl = get(ax(i), 'XLabel'); set(xl, 'Color', 'k');
        yl = get(ax(i), 'YLabel'); set(yl, 'Color', 'k');
    end
    lgd = findall(gcf, 'Type', 'legend');
    set(lgd, 'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'k');
end
