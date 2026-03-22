%% Risk-Aware Power System - Algorithm 2 Multi-Case Study
% This script evaluates Algorithm 2 for multiple heterogeneous initial
% conditions beta_i'(0):
%   1) consensus estimate of sum(beta_i')
%   2) distributed update of beta_i' toward (1-n)/4
%   3) one-shot bid alpha_i* from current beta_i'
%
% Figures are saved to ./figures using names referenced by main.tex.

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

n  = 3;
a2 = 0.1;
a1 = 1;
a0 = 0; %#ok<NASGU>

x0  = [120; 100; 80];
xL0 = 400;

sigma_r2 = 0;
sigma_L2 = 0;
sigma_rL = 0;

% Operator-only communication graph (connected)
Adj_ops = ones(n) - eye(n);
L_ops   = diag(sum(Adj_ops, 2)) - Adj_ops;

% Consensus + outer-loop parameters
epsilon_c     = 0.2;      % consensus gain (< 1/max_degree)
max_cons_iter = 150;
cons_tol      = 1e-10;

epsilon_beta  = 0.25;     % require 0 < n*epsilon_beta < 2 for linear model
delta_beta    = 1e-4;     % keeps beta_i' > -1/4
beta_tol      = 1e-8;
f_tol         = 1e-8;
max_outer_iter = 300;

% If true: enforce no-overbidding via alpha_i in [0,1]
enforce_alpha_clamp = false;

save_figs  = true;
plot_figs  = true;   % set false to skip all plotting
close_figs = true;   % set true to close all figures after saving
fig_out_dir = 'figures';
plot_tail_buffer = 12;     % extra points after first practical convergence
plot_close_tol = 1e-4;     % closeness to final value for practical convergence

%% ========================================================================
%                         TEST CASES (beta_i'(0))
% =========================================================================

test_cases = {
    'Case I',   [0; 0; 0],            'Risk-neutral';
    'Case II',  [-1; -1; -1],         'Risk-seeking';
    'Case III', [1; 1; 1],            'Risk-averse';
    'Case IV',  [-1; 1; 1],           'Mixed';
    'Case V',   [-0.25; 0; 0.20],     'Mixed (mild)';
    'Case VI',  [-0.25; 0; 0],        'Mixed (mild)';
};

num_cases = size(test_cases, 1);
case_names = test_cases(:, 1);
case_betas = test_cases(:, 2);
case_profiles = test_cases(:, 3);

%% ========================================================================
%                         DERIVED QUANTITIES              5
% =========================================================================

xr0 = sum(x0);
eta_true = x0 / xr0;
B_target = (1 - n) / 4;
g_max = 0.25 * xr0^2;

% Distributed estimate of xr0 and eta (operator-only, for this test script)
z_xr = x0;
for t = 1:max_cons_iter
    z_new = z_xr - epsilon_c * L_ops * z_xr;
    if norm(z_new - z_xr) < cons_tol
        z_xr = z_new;
        break;
    end
    z_xr = z_new;
end
xr0_est = n * mean(z_xr);
eta_est = x0 / xr0_est;

% Centralized references
[alpha_RN, f_RN, Pi_RN] = nash_eq(zeros(n,1), n, a2, a1, x0, xL0, xr0, eta_true, sigma_r2, sigma_L2, sigma_rL);
beta_prime_sym = (1 - n) / (4 * n);
[alpha_sym, f_sym, Pi_sym] = nash_eq(beta_prime_sym * ones(n,1), n, a2, a1, x0, xL0, xr0, eta_true, sigma_r2, sigma_L2, sigma_rL);

fprintf('==========================================================\n');
fprintf('   ALGORITHM 2 MULTI-CASE STUDY (NON-EQUAL beta_i''(0))\n');
fprintf('==========================================================\n\n');
fprintf('n=%d, a2=%.2f, x0=[%s], xr0=%.1f\n', n, a2, num2str(x0'), xr0);
fprintf('eta_true=[%s], eta_est=[%s]\n', num2str(eta_true', '%.6f '), num2str(eta_est', '%.6f '));
fprintf('Target sum(beta'')=(1-n)/4=%.6f\n\n', B_target);

%% ========================================================================
%                         RUN ALL CASES
% =========================================================================

all_results = cell(num_cases, 1);

for tc = 1:num_cases
    beta0 = case_betas{tc};
    beta_k = max(beta0, -0.25 + delta_beta);

    hist.beta = zeros(n, max_outer_iter);
    hist.Bsum = zeros(1, max_outer_iter);
    hist.eB = zeros(1, max_outer_iter);
    hist.alpha_unc = zeros(n, max_outer_iter);
    hist.alpha = zeros(n, max_outer_iter);
    hist.f = zeros(1, max_outer_iter);
    hist.S = zeros(1, max_outer_iter);
    hist.g = zeros(1, max_outer_iter);
    hist.Pi = zeros(n, max_outer_iter);
    hist.Pi_clp = zeros(n, max_outer_iter);
    hist.sat_count = zeros(1, max_outer_iter);
    hist.sat_per_op = zeros(n, max_outer_iter);

    converged = false;
    iters = max_outer_iter;

    for k = 1:max_outer_iter
        hist.beta(:,k) = beta_k;

        % Consensus estimate of aggregate risk level B(k)=sum beta_i'(k)
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
        eB = B_est - B_target;

        % One-shot bidding induced by current beta
        alpha_unc = 0.5 * (1 + 4*beta_k) ./ eta_est;
        if enforce_alpha_clamp
            alpha_k = max(0, min(1, alpha_unc));
        else
            alpha_k = alpha_unc;
        end

        f_k = sum(alpha_k .* eta_est);
        S_k = f_k - 0.5;
        g_k = f_k * (1 - f_k) * xr0^2;
        Pi_k = profit(alpha_k, a2, a1, x0, xL0, xr0, eta_est, sigma_r2, sigma_L2, sigma_rL);
        alpha_clp_k = max(0, min(1, alpha_unc));
        Pi_clp_k = profit(alpha_clp_k, a2, a1, x0, xL0, xr0, eta_est, sigma_r2, sigma_L2, sigma_rL);

        hist.Bsum(k) = B_est;
        hist.eB(k) = eB;
        hist.alpha_unc(:,k) = alpha_unc;
        hist.alpha(:,k) = alpha_k;
        hist.f(k) = f_k;
        hist.S(k) = S_k;
        hist.g(k) = g_k;
        hist.Pi(:,k) = Pi_k;
        hist.Pi_clp(:,k) = Pi_clp_k;
        hist.sat_count(k) = sum(alpha_unc > 1);
        hist.sat_per_op(:,k) = alpha_unc > 1;

        if abs(eB) < beta_tol && abs(S_k) < f_tol
            converged = true;
            iters = k;
            break;
        end

        % Algorithm 2 update on beta only (aggregate-error correction)
        beta_k = max(beta_k - epsilon_beta * eB, -0.25 + delta_beta);
    end

    hist.beta = hist.beta(:,1:iters);
    hist.Bsum = hist.Bsum(1:iters);
    hist.eB = hist.eB(1:iters);
    hist.alpha_unc = hist.alpha_unc(:,1:iters);
    hist.alpha = hist.alpha(:,1:iters);
    hist.f = hist.f(1:iters);
    hist.S = hist.S(1:iters);
    hist.g = hist.g(1:iters);
    hist.Pi = hist.Pi(:,1:iters);
    hist.Pi_clp = hist.Pi_clp(:,1:iters);
    hist.sat_count = hist.sat_count(1:iters);
    hist.sat_per_op = hist.sat_per_op(:,1:iters);
    hist.beta0 = beta0;
    hist.beta_final = beta_k;
    hist.converged = converged;
    hist.iters = iters;
    idx_plot = find(abs(hist.f - hist.f(end)) < plot_close_tol & ...
                    abs(hist.Bsum - hist.Bsum(end)) < plot_close_tol, 1, 'first');
    if isempty(idx_plot)
        hist.k_plot = iters;
    else
        hist.k_plot = min(iters, idx_plot + plot_tail_buffer);
    end

    all_results{tc} = hist;

    fprintf('%-7s  iters=%3d  sum(beta'')=% .6f  f=% .6f  sum(Pi)=%.2f  sat=%d\n', ...
        case_names{tc}, iters, sum(hist.beta_final), hist.f(end), sum(hist.Pi(:,end)), hist.sat_count(end));
end

%% ========================================================================
%                            SUMMARY TABLE
% =========================================================================

fprintf('\n----------------------------------------------------------\n');
fprintf('References: RN sum(Pi)=%.2f, Pareto-central sum(Pi)=%.2f\n', sum(Pi_RN), sum(Pi_sym));
fprintf('----------------------------------------------------------\n');
fprintf('%-8s | %5s | %11s | %9s | %10s | %10s\n', 'Case', 'Iters', 'sum(beta'')', 'f(final)', 'g(final)', 'sum(Pi)');
fprintf('%s\n', repmat('-',1,72));
for tc = 1:num_cases
    r = all_results{tc};
    fprintf('%-8s | %5d | %11.6f | %9.6f | %10.2f | %10.2f\n', ...
        case_names{tc}, r.iters, sum(r.beta_final), r.f(end), r.g(end), sum(r.Pi(:,end)));
end

%% ========================================================================
%                           VISUALIZATION
% =========================================================================

if ~plot_figs
    fprintf('\nplot_figs=false: skipping all figures.\n');
    return;
end

if save_figs && ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

colors = lines(num_cases);
k_plot_global = max(cellfun(@(r) r.k_plot, all_results));

% 1) sum beta convergence
fig1 = figure('Position', [50 50 900 520]); hold on;
for tc = 1:num_cases
    r = all_results{tc};
    kshow = r.k_plot;
    plot(1:kshow, r.Bsum(1:kshow), 'LineWidth', 2, 'Color', colors(tc,:), 'DisplayName', case_names{tc});
end
yline(B_target, 'k--', 'LineWidth', 1.8, 'DisplayName', 'target (1-n)/4');
xlabel('Iteration k'); ylabel('\Sigma\beta_i''(k)');
title('Convergence of \Sigma\beta_i'' for all test cases');
xlim([1, max(2, k_plot_global)]);
legend('Location', 'best'); grid on;
style_figure(fig1);
if save_figs, exportgraphics(fig1, fullfile(fig_out_dir, 'sum_beta_convergence_all_cases.png'), 'Resolution', 150); end

% 3) alpha individual
fig3 = figure('Position', [100 70 1250 820]);
for tc = 1:num_cases
    subplot(2,3,tc);
    r = all_results{tc};
    kshow = r.k_plot;
    plot(1:kshow, r.alpha(:,1:kshow)', 'LineWidth', 1.4); hold on;
    yline(1, 'k--', 'LineWidth', 1.0);
    title(sprintf('%s (%s)', case_names{tc}, case_profiles{tc}));
    xlim([1, max(2, kshow)]);
    xlabel('k'); ylabel('\alpha_i'); grid on;
end
sg3 = sgtitle('Individual \alpha_i convergence for all test cases');
set(sg3, 'Color', 'k');
style_figure(fig3);
if save_figs, exportgraphics(fig3, fullfile(fig_out_dir, 'alpha_individual_all_cases.png'), 'Resolution', 150); end

% 4) f convergence
fig4 = figure('Position', [120 80 900 520]); hold on;
for tc = 1:num_cases
    r = all_results{tc};
    kshow = r.k_plot;
    plot(1:kshow, r.f(1:kshow), 'LineWidth', 2, 'Color', colors(tc,:), 'DisplayName', case_names{tc});
end
yline(0.5, 'k--', 'LineWidth', 1.8, 'DisplayName', 'target 0.5');
xlabel('Iteration k'); ylabel('f(\alpha)');
title('Convergence of f(\alpha) to Pareto target 0.5');
xlim([1, max(2, k_plot_global)]);
legend('Location', 'best'); grid on;
style_figure(fig4);
if save_figs, exportgraphics(fig4, fullfile(fig_out_dir, 'f_convergence_all_cases.png'), 'Resolution', 150); end

% 5) S convergence (log)
fig5 = figure('Position', [140 90 900 520]); hold on;
for tc = 1:num_cases
    r = all_results{tc};
    kshow = r.k_plot;
    semilogy(1:kshow, max(abs(r.S(1:kshow)), 1e-12), 'LineWidth', 2, 'Color', colors(tc,:), 'DisplayName', case_names{tc});
end
xlabel('Iteration k'); ylabel('|S|=|f-0.5|');
title('Convergence of error signal |S| (log scale)');
xlim([1, max(2, k_plot_global)]);
legend('Location', 'best'); grid on;
style_figure(fig5);
if save_figs, exportgraphics(fig5, fullfile(fig_out_dir, 'S_convergence_all_cases.png'), 'Resolution', 150); end

% 6) g convergence
fig6 = figure('Position', [160 100 900 520]); hold on;
for tc = 1:num_cases
    r = all_results{tc};
    kshow = r.k_plot;
    plot(1:kshow, r.g(1:kshow), 'LineWidth', 2, 'Color', colors(tc,:), 'DisplayName', case_names{tc});
end
yline(g_max, 'k--', 'LineWidth', 1.8, 'DisplayName', 'g_{max}');
xlabel('Iteration k'); ylabel('g(\alpha)=f(1-f)(x_r^o)^2');
title('Convergence of market efficiency g(\alpha)');
xlim([1, max(2, k_plot_global)]);
legend('Location', 'best'); grid on;
style_figure(fig6);
if save_figs, exportgraphics(fig6, fullfile(fig_out_dir, 'g_convergence_all_cases.png'), 'Resolution', 150); end

% 7) aggregate profit convergence
fig7 = figure('Position', [180 110 900 520]); hold on;
for tc = 1:num_cases
    r = all_results{tc};
    kshow = r.k_plot;
    plot(1:kshow, sum(r.Pi(:,1:kshow),1), 'LineWidth', 2, 'Color', colors(tc,:), 'DisplayName', case_names{tc});
end
yline(sum(Pi_sym), 'k--', 'LineWidth', 1.8, 'DisplayName', 'Pareto central');
xlabel('Iteration k'); ylabel('\Sigma\Pi_i');
title('Convergence of aggregate profit');
xlim([1, max(2, k_plot_global)]);
legend('Location', 'best'); grid on;
style_figure(fig7);
if save_figs, exportgraphics(fig7, fullfile(fig_out_dir, 'profit_convergence_all_cases.png'), 'Resolution', 150); end

% 8) individual profit convergence
fig8 = figure('Position', [200 120 1250 820]);
for tc = 1:num_cases
    subplot(2,3,tc);
    r = all_results{tc};
    kshow = r.k_plot;
    plot(1:kshow, r.Pi(:,1:kshow)', 'LineWidth', 1.4);
    title(sprintf('%s (%s)', case_names{tc}, case_profiles{tc}));
    xlim([1, max(2, kshow)]);
    xlabel('k'); ylabel('\Pi_i'); grid on;
end
sg8 = sgtitle('Individual profit convergence for all test cases');
set(sg8, 'Color', 'k');
style_figure(fig8);
if save_figs, exportgraphics(fig8, fullfile(fig_out_dir, 'profit_individual_all_cases.png'), 'Resolution', 150); end

% 9) Overbidding bar chart: final alpha values and # overbidding operators per case
% Show which cases exhibit overbidding at convergence and what the bids look like.
op_colors_local = lines(n);
case_cats = categorical(case_names, case_names);

% Gather final alpha (unclamped) and per-operator overbidding flag for each case
alpha_final_all  = zeros(n, num_cases);   % unclamped final alpha per case
sat_final_all    = zeros(n, num_cases);   % 1 = overbidding at final iter
for tc = 1:num_cases
    r = all_results{tc};
    alpha_final_all(:, tc)  = r.alpha_unc(:, end);
    sat_final_all(:, tc)    = double(r.alpha_unc(:, end) > 1);
end

fig9 = figure('Position', [220 130 720 460]);
b_af = bar(case_cats, alpha_final_all', 'grouped');
for i = 1:n
    b_af(i).FaceColor = op_colors_local(i,:);
    b_af(i).DisplayName = sprintf('Op%d', i);
end
hold on;
hl9 = yline(1, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Capacity limit');
xlabel('Test Case'); ylabel('\alpha_i (final, unclamped)');
title('Final Bids at Convergence (Unclamped)');
legend([b_af, hl9], 'Location', 'best', 'FontSize', 7);
grid on;
style_figure(fig9);
if save_figs, exportgraphics(fig9, fullfile(fig_out_dir, 'mc_overbid_count.png'), 'Resolution', 150); end

if save_figs
    fprintf('\nSaved figures to %s/\n', fig_out_dir);
end

if close_figs
    close all;
    fprintf('All figures closed (close_figs=true).\n');
end

%% ========================================================================
%                         HELPER FUNCTIONS
% =========================================================================

function [alpha_star, f, Pi] = nash_eq(beta_prime, n, a2, a1, x0, xL0, xr0, eta, sigma_r2, sigma_L2, sigma_rL)
    B_sum = sum(beta_prime);
    f = (n + 4*B_sum) / (1 + n + 4*B_sum);
    alpha_star = (1 + 4*beta_prime) .* (1-f) ./ eta;
    alpha_star = max(0, min(1, alpha_star));
    f = sum(alpha_star .* eta);
    Pi = profit(alpha_star, a2, a1, x0, xL0, xr0, eta, sigma_r2, sigma_L2, sigma_rL);
end

function Pi = profit(alpha, a2, a1, x0, xL0, xr0, eta, sigma_r2, sigma_L2, sigma_rL)
    f  = sum(alpha .* eta);
    Pi = 2*a2 * alpha .* x0 * xr0 .* (1-f) ...
       + a1 * x0 ...
       + 2*a2 * x0 .* (xL0 - xr0);
    if sigma_rL ~= 0 || sigma_r2 > 0
        Pi = Pi + 2*a2*(sigma_rL - sigma_r2/length(x0)) * ones(length(x0),1);
    end
end

function style_figure(fig_h)
    figure(fig_h); %#ok<LFIG>
    set(fig_h, 'Color', 'w', 'InvertHardcopy', 'off');
    ax = findall(fig_h, 'Type', 'axes');
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    for i = 1:numel(ax)
        t = get(ax(i), 'Title');
        set(t, 'Color', 'k');
        xl = get(ax(i), 'XLabel');
        yl = get(ax(i), 'YLabel');
        set(xl, 'Color', 'k');
        set(yl, 'Color', 'k');
    end
    lgd = findall(fig_h, 'Type', 'legend');
    set(lgd, 'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'k');
    tx = findall(fig_h, 'Type', 'text');
    set(tx, 'Color', 'k');
    st = findall(fig_h, 'Tag', 'sgtitle');
    set(st, 'Color', 'k');
end
