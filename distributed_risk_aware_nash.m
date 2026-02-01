%% Distributed Risk-Aware Nash Optimization Algorithm
% Implementation of Algorithm 1 from the paper
% "Distributed Risk-Aware Bidding Strategy for Incorporating Renewable
% Generation into Real-Time Electricity Market"
%
% SINGLE TEST CASE VERSION - for debugging and parameter exploration
% For multiple test cases comparison, use risk_aware_power_system.m
%
% Two modes:
%   1. Centralized: Ground truth solution using global information
%   2. Distributed: Consensus-based algorithm using only local communication

clc; clear; close all;

%% ========================================================================
%                           USER PARAMETERS
% =========================================================================

% Number of renewable operators
n = 3;

% Cost coefficients (utility company)
a2 = 0.1;
a1 = 1;
a0 = 0;

% Forecasted renewable generations x_i^o (can be asymmetric)
x0 = [120; 100; 80];  % MW

% Forecasted load x_{n+1}^o
xL0 = 400;  % MW

% Covariance parameters (set to zero for deterministic analysis)
sigma_r2 = 0;      % variance of aggregate renewable
sigma_L2 = 0;      % variance of load
sigma_rL = 0;      % covariance between renewable and load

% Initial bidding fractions (user choice)
alpha_init = [0.5; 0.5; 0.5];

% Initial transformed risk parameters beta_i' = a2 * beta_i (user choice)
% Negative = risk-seeking, Zero = risk-neutral, Positive = risk-averse
beta_prime_init = [-0.2; -0.23; -0.25];  % Mixed case

% Communication graph (adjacency matrix) - must be connected
% Ring topology: 1 -- 2 -- 3 -- 1
Adj = [0 1 1;
       1 0 1;
       1 1 0];  % Fully connected

% Algorithm parameters
epsilon1 = 0.2;      % Consensus step size (< 1/max_degree)
epsilon2 = 0.05;     % Beta update step size (increased for faster convergence)
tau = 1e-6;          % Stopping tolerance for |S| = |f - 0.5|
alpha_tol = 1e-6;    % Stopping tolerance for alpha convergence
min_iter = 10;       % Minimum iterations before checking convergence
max_outer_iter = 1000; % Max iterations for outer loop (increased)
max_consensus_iter = 100;  % Max iterations for consensus
consensus_tol = 1e-10;  % Consensus convergence tolerance

% Visualization parameters
save_figs = false;           % Set to true to save figures to figs/ folder
use_subplots = true;        % Set to true for one big subplot, false for individual figures

%% ========================================================================
%                         DERIVED QUANTITIES
% =========================================================================

xr0 = sum(x0);           % Aggregate renewable forecast
eta = x0 / xr0;          % Forecast ratios eta_i = x_i^o / x_r^o
gamma = x0' ./ x0;       % gamma_{ij} = x_j^o / x_i^o (n x n matrix)

% Laplacian matrix for consensus
Deg = diag(sum(Adj, 2));
L = Deg - Adj;

fprintf('=========================================================\n');
fprintf('   DISTRIBUTED RISK-AWARE NASH OPTIMIZATION\n');
fprintf('             SINGLE TEST CASE (DEBUG MODE)\n');
fprintf('=========================================================\n\n');
fprintf('Problem Parameters:\n');
fprintf('  n = %d operators\n', n);
fprintf('  a2 = %.2f, a1 = %.2f, a0 = %.2f\n', a2, a1, a0);
fprintf('  x0 = [%s] MW\n', num2str(x0'));
fprintf('  xL0 = %.2f MW\n', xL0);
fprintf('  x_r^o = %.2f MW\n', xr0);
fprintf('  eta = [%s]\n\n', num2str(eta', '%.4f '));

% Warn if initial beta' is at or below stability boundary
if any(beta_prime_init <= -0.25)
    warning('Initial beta_prime contains values at or below -0.25 (stability boundary).');
    fprintf('  Values will be clamped to > -0.25 during iteration.\n\n');
end

% Target sum of beta' for Pareto optimal f = 0.5
beta_prime_sum_target = (1 - n) / 4;
fprintf('Target sum(beta'') for f=0.5: %.4f\n', beta_prime_sum_target);
fprintf('Initial sum(beta''): %.4f\n\n', sum(beta_prime_init));

%% ========================================================================
%                    PART 1: CENTRALIZED SOLUTION (GROUND TRUTH)
% =========================================================================

fprintf('=========================================================\n');
fprintf('              PART 1: CENTRALIZED SOLUTION\n');
fprintf('=========================================================\n\n');

% --- 1a. Risk-Neutral Nash Equilibrium (beta_i = 0) ---
fprintf('--- 1a. Risk-Neutral Nash Equilibrium (beta_i = 0) ---\n');
[alpha_RN, f_RN, Pi_RN, J_RN] = compute_nash_centralized(zeros(n,1), ...
    n, a2, x0, xL0, xr0, eta, gamma, sigma_r2, sigma_L2, sigma_rL);

fprintf('alpha_RN = [%s]\n', num2str(alpha_RN', '%.6f '));
fprintf('f_RN = %.6f (Pareto optimal: 0.5)\n', f_RN);
fprintf('Aggregate profit Pi_RN = %.4f\n', sum(Pi_RN));
fprintf('g_RN = f*(1-f)*x_r^o^2 = %.4f (max possible: %.4f)\n\n', ...
    f_RN*(1-f_RN)*xr0^2, 0.25*xr0^2);

% --- 1b. Symmetric Pareto-Optimal Strategy (Corollary 5) ---
fprintf('--- 1b. Symmetric Pareto-Optimal Strategy ---\n');
beta_prime_sym = (1 - n) / (4 * n);  % Eq. (4.6)
beta_sym = beta_prime_sym / a2;      % Convert back to beta

fprintf('Theoretical beta_sym'' = (1-n)/(4n) = %.6f\n', beta_prime_sym);
fprintf('Theoretical beta_sym = beta_sym''/a2 = %.6f\n', beta_sym);

[alpha_sym, f_sym, Pi_sym, J_sym] = compute_nash_centralized(...
    beta_prime_sym * ones(n,1), n, a2, x0, xL0, xr0, eta, gamma, ...
    sigma_r2, sigma_L2, sigma_rL);

fprintf('alpha_sym = [%s]\n', num2str(alpha_sym', '%.6f '));
fprintf('f_sym = %.6f (target: 0.5)\n', f_sym);
fprintf('Aggregate profit Pi_sym = %.4f\n', sum(Pi_sym));
fprintf('g_sym = f*(1-f)*x_r^o^2 = %.4f\n\n', f_sym*(1-f_sym)*xr0^2);

% --- 1c. User-specified initial parameters ---
fprintf('--- 1c. Nash with User Initial Parameters ---\n');
[alpha_opt, f_user, Pi_user, J_user] = compute_nash_centralized(...
    beta_prime_init, n, a2, x0, xL0, xr0, eta, gamma, ...
    sigma_r2, sigma_L2, sigma_rL);

fprintf('Initial beta'' = [%s]\n', num2str(beta_prime_init', '%.6f '));
fprintf('alpha_opt = [%s]\n', num2str(alpha_opt', '%.6f '));
fprintf('f_user = %.6f\n', f_user);
fprintf('Aggregate profit = %.4f\n\n', sum(Pi_user));

%% ========================================================================
%                    PART 2: DISTRIBUTED ALGORITHM
% =========================================================================

fprintf('=========================================================\n');
fprintf('              PART 2: DISTRIBUTED ALGORITHM\n');
fprintf('=========================================================\n\n');

% Initialize variables
alpha = alpha_init;
beta_prime = beta_prime_init;

% History for plotting
history.alpha = zeros(n, max_outer_iter);
history.beta_prime = zeros(n, max_outer_iter);
history.f = zeros(1, max_outer_iter);
history.S = zeros(1, max_outer_iter);
history.Pi = zeros(n, max_outer_iter);
history.J = zeros(n, max_outer_iter);
history.alpha_change = zeros(1, max_outer_iter);

% --- Step 1: One-time consensus for aggregate quantities ---
fprintf('Step 1: One-time consensus for n, x_r^o, eta_i\n');

% Compute n via consensus (Eq. 5.3)
z = zeros(n, 1); z(1) = 1;
n_est = 1 / run_consensus(z, L, epsilon1, max_consensus_iter, consensus_tol);
fprintf('  Estimated n = %.4f (true = %d)\n', n_est, n);

% Compute x_r^o via consensus (Eq. 5.4)
z = x0;
z_bar = run_consensus(z, L, epsilon1, max_consensus_iter, consensus_tol);
xr0_est = n_est * z_bar;
fprintf('  Estimated x_r^o = %.4f (true = %.2f)\n', xr0_est, xr0);

% Compute eta_i locally (Eq. 5.5)
eta_est = x0 / xr0_est;
fprintf('  Estimated eta = [%s]\n\n', num2str(eta_est', '%.4f '));

% --- Step 2: Iterative optimization ---
fprintf('Step 2: Iterative optimization loop\n');
fprintf('  Initial alpha = [%s]\n', num2str(alpha', '%.4f '));
fprintf('  Initial beta'' = [%s]\n\n', num2str(beta_prime', '%.4f '));

converged = false;
alpha_prev = alpha;  % Track previous alpha for convergence check

for k = 1:max_outer_iter
    % Store history
    history.alpha(:, k) = alpha;
    history.beta_prime(:, k) = beta_prime;

    % Compute f via consensus (Eq. 5.6)
    z = alpha .* x0;
    z_bar = run_consensus(z, L, epsilon1, max_consensus_iter, consensus_tol);
    f = n_est * z_bar / xr0_est;
    history.f(k) = f;

    % Compute error signal S (Eq. 5.8)
    S = f - 0.5;
    history.S(k) = S;

    % Compute profits and objectives at current state
    [~, ~, Pi_k, J_k] = compute_nash_centralized(beta_prime, n, a2, ...
        x0, xL0, xr0, eta, gamma, sigma_r2, sigma_L2, sigma_rL);
    % Use current alpha for actual profit calculation
    Pi_actual = compute_profit(alpha, n, a2, a1, x0, xL0, xr0, eta, ...
        sigma_r2, sigma_L2, sigma_rL);
    J_actual = Pi_actual - beta_prime .* 4*a2*(xr0^2)*(1-f)^2;
    history.Pi(:, k) = Pi_actual;
    history.J(:, k) = J_actual;

    % Update alpha using best response (derived from Nash condition)
    % alpha_i = (1 + 4*beta_i') * (1 - f) / eta_i
    % But clamp to [0, 1] and use damping for stability
    alpha_target = (1 + 4*beta_prime) .* (1 - f) ./ eta_est;
    alpha_target = max(0, min(1, alpha_target));  % Clamp to [0, 1]
    alpha_new = 0.5 * alpha + 0.5 * alpha_target;  % Damped update

    % Check stopping criterion (Eq. 5.10)
    % Only check after minimum iterations AND require both:
    % 1. |S| < tau (f is at Pareto target)
    % 2. alpha has stabilized (at Nash equilibrium for current beta)
    alpha_change = norm(alpha_new - alpha);
    history.alpha_change(k) = alpha_change;

    if k >= min_iter && abs(S) < tau && alpha_change < alpha_tol
        fprintf('  Converged at iteration %d: |S| = %.2e, |delta_alpha| = %.2e\n', ...
            k, abs(S), alpha_change);
        converged = true;
        break;
    end

    % Apply alpha update
    alpha_prev = alpha;
    alpha = alpha_new;

    % Update beta' using gradient descent (Eq. 5.9)
    beta_prime = beta_prime - epsilon2 * S;

    % Ensure beta' > -1/4 for stability (Theorem 1)
    % Use -0.2499 to allow starting at -0.25 without immediate clamping issues
    beta_prime = max(beta_prime, -0.2499);

    % Progress display
    if mod(k, 50) == 0
        fprintf('  Iter %d: f = %.6f, S = %.6f, |d_alpha| = %.2e, sum(Pi) = %.4f\n', ...
            k, f, S, alpha_change, sum(Pi_actual));
    end
end

if ~converged
    fprintf('  Did not converge within %d iterations\n', max_outer_iter);
end

% Trim history
history.alpha = history.alpha(:, 1:k);
history.beta_prime = history.beta_prime(:, 1:k);
history.f = history.f(1:k);
history.S = history.S(1:k);
history.Pi = history.Pi(:, 1:k);
history.J = history.J(:, 1:k);
history.alpha_change = history.alpha_change(1:k);

%% ========================================================================
%                         FINAL RESULTS
% =========================================================================

fprintf('\n=========================================================\n');
fprintf('                    FINAL RESULTS\n');
fprintf('=========================================================\n\n');

% Final values from distributed algorithm
alpha_final = history.alpha(:, end);
beta_prime_final = history.beta_prime(:, end);
f_final = history.f(end);
Pi_final = history.Pi(:, end);
J_final = history.J(:, end);

fprintf('--- Distributed Algorithm Output ---\n');
fprintf('Pareto-optimal beta_i'' = [%s]\n', num2str(beta_prime_final', '%.6f '));
fprintf('Pareto-optimal beta_i  = [%s]\n', num2str((beta_prime_final/a2)', '%.6f '));
fprintf('Pareto-optimal alpha_i = [%s]\n', num2str(alpha_final', '%.6f '));
fprintf('Final f = %.6f (target: 0.5)\n', f_final);
fprintf('Final |S| = %.2e\n\n', abs(history.S(end)));

fprintf('--- Profit and Objective ---\n');
fprintf('Individual profits Pi_i = [%s]\n', num2str(Pi_final', '%.4f '));
fprintf('Aggregate profit sum(Pi) = %.4f\n', sum(Pi_final));
fprintf('Individual objectives J_i = [%s]\n', num2str(J_final', '%.4f '));
fprintf('Aggregate objective sum(J) = %.4f\n\n', sum(J_final));

fprintf('--- Comparison with Centralized ---\n');
fprintf('Centralized symmetric beta'' = %.6f\n', beta_prime_sym);
fprintf('Distributed mean beta''     = %.6f\n', mean(beta_prime_final));
fprintf('Centralized symmetric f     = %.6f\n', f_sym);
fprintf('Distributed f               = %.6f\n', f_final);
fprintf('Centralized symmetric g     = %.4f\n', f_sym*(1-f_sym)*xr0^2);
fprintf('Distributed g               = %.4f\n', f_final*(1-f_final)*xr0^2);
fprintf('Maximum possible g          = %.4f\n\n', 0.25*xr0^2);

fprintf('--- Sum(beta'') Convergence ---\n');
beta_prime_sum_target = (1 - n) / 4;
fprintf('Target sum(beta'') for f=0.5: %.6f\n', beta_prime_sum_target);
fprintf('Initial sum(beta''):         %.6f\n', sum(beta_prime_init));
fprintf('Final sum(beta''):           %.6f\n', sum(beta_prime_final));
fprintf('Distance to target:          %.6f\n\n', abs(sum(beta_prime_final) - beta_prime_sum_target));

%% ========================================================================
%                         VISUALIZATION
% =========================================================================

% Create figs folder if saving is enabled
if save_figs && ~exist('figs', 'dir')
    mkdir('figs');
end

% Generate configuration string for filenames
config_str = sprintf('n%d_alpha%.1f-%.1f-%.1f_beta%.2f-%.2f-%.2f', ...
    n, alpha_init(1), alpha_init(2), alpha_init(3), ...
    beta_prime_init(1), beta_prime_init(2), beta_prime_init(3));

% Compute optimal sum(beta') for Pareto condition
% For general case, sum(beta') should satisfy the Pareto condition
% For symmetric case: n * beta_sym' = n * (1-n)/(4n) = (1-n)/4
beta_prime_sum_opt = n * beta_prime_sym;  % Optimal sum of beta'

% Compute sum(beta') history
beta_prime_sum_history = sum(history.beta_prime, 1);  % Convert back to beta'

% Common figure settings
fig_width = 600;
fig_height = 450;

if use_subplots
    % === SUBPLOT MODE: One big figure with all plots ===
    fig_main = figure('Position', [50, 50, 1400, 900], 'Name', 'Distributed Risk-Aware Nash');

    % --- Plot 1: Alpha convergence ---
    subplot(2, 3, 1);
    plot(1:k, history.alpha', 'LineWidth', 1.5);
    hold on;
    yline(0.5, 'k--', 'LineWidth', 1.5);
    xlabel('Iteration k');
    ylabel('\alpha_i(k)');
    title('Bidding Fraction Convergence');
    legend_str = cell(n+1, 1);
    for i = 1:n
        legend_str{i} = sprintf('Op %d', i);
    end
    legend_str{n+1} = 'Target';
    legend(legend_str, 'Location', 'best', 'FontSize', 8);
    grid on;

    % --- Plot 2: Beta' convergence with sum ---
    subplot(2, 3, 2);
    plot(1:k, history.beta_prime', 'LineWidth', 1.5);
    hold on;
    plot(1:k, beta_prime_sum_history, 'k-', 'LineWidth', 2);
    yline(beta_prime_sym, 'r--', 'LineWidth', 1.5);
    yline(beta_prime_sum_opt, 'k--', 'LineWidth', 1.5);
    xlabel('Iteration k');
    ylabel('\beta_i''(k)');
    title('Risk Parameter Convergence');
    legend_str = cell(n+3, 1);
    for i = 1:n
        legend_str{i} = sprintf('Op %d', i);
    end
    legend_str{n+1} = '\Sigma\beta_i''';
    legend_str{n+2} = sprintf('\\beta_{sym}^* (%.3f)', beta_prime_sym);
    legend_str{n+3} = sprintf('\\Sigma\\beta^* (%.3f)', beta_prime_sum_opt);
    legend(legend_str, 'Location', 'best', 'FontSize', 8);
    grid on;

    % --- Plot 3: f and S convergence ---
    subplot(2, 3, 3);
    yyaxis left;
    plot(1:k, history.f, 'b-', 'LineWidth', 1.5);
    hold on;
    yline(0.5, 'b--', 'LineWidth', 1);
    ylabel('f(k)');
    yyaxis right;
    semilogy(1:k, abs(history.S), 'r-', 'LineWidth', 1.5);
    hold on;
    semilogy(1:k, history.alpha_change, 'm-', 'LineWidth', 1.5);
    yline(tau, 'r--', 'LineWidth', 1);
    ylabel('|S|, |\Delta\alpha|');
    xlabel('Iteration k');
    title('Convergence Metrics');
    legend('f', 'f=0.5', '|S|', '|\Delta\alpha|', '\tau', 'Location', 'best', 'FontSize', 8);
    grid on;

    % --- Plot 4: Individual profits ---
    subplot(2, 3, 4);
    plot(1:k, history.Pi', 'LineWidth', 1.5);
    xlabel('Iteration k');
    ylabel('\Pi_i(k)');
    title('Individual Profit Convergence');
    legend_str = cell(n, 1);
    for i = 1:n
        legend_str{i} = sprintf('Op %d', i);
    end
    legend(legend_str, 'Location', 'best', 'FontSize', 8);
    grid on;

    % --- Plot 5: Individual objectives ---
    subplot(2, 3, 5);
    plot(1:k, history.J', 'LineWidth', 1.5);
    xlabel('Iteration k');
    ylabel('J_i(k)');
    title('Individual Objective Convergence');
    legend(legend_str, 'Location', 'best', 'FontSize', 8);
    grid on;

    % --- Plot 6: Aggregate metrics ---
    subplot(2, 3, 6);
    g_history = history.f .* (1 - history.f) * xr0^2;
    plot(1:k, sum(history.Pi, 1), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(1:k, g_history, 'r-', 'LineWidth', 1.5);
    yline(0.25*xr0^2, 'r--', 'LineWidth', 1);
    yline(sum(Pi_sym), 'b--', 'LineWidth', 1);
    xlabel('Iteration k');
    ylabel('Value');
    title('Aggregate Metrics Convergence');
    legend('\Sigma\Pi_i', 'g', 'g_{max}', '\Pi_{sym}', 'Location', 'best', 'FontSize', 8);
    grid on;

    sgtitle('Distributed Risk-Aware Nash Optimization', 'FontSize', 14, 'FontWeight', 'bold');

    if save_figs
        saveas(fig_main, fullfile('figs', ['all_plots_' config_str '.png']));
        fprintf('Figure saved to figs/all_plots_%s.png\n', config_str);
    end

else
    % === INDIVIDUAL FIGURE MODE ===

    % --- Plot 1: Alpha convergence ---
    fig1 = figure('Position', [50, 50, fig_width, fig_height]);
    plot(1:k, history.alpha', 'LineWidth', 1.5);
    hold on;
    yline(0.5, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Pareto target');
    xlabel('Iteration k');
    ylabel('\alpha_i(k)');
    title('Bidding Fraction Convergence');
    legend_str = cell(n+1, 1);
    for i = 1:n
        legend_str{i} = sprintf('Operator %d', i);
    end
    legend_str{n+1} = 'Pareto target (\alpha=0.5)';
    legend(legend_str, 'Location', 'best');
    grid on;
    if save_figs
        saveas(fig1, fullfile('figs', ['alpha_convergence_' config_str '.png']));
    end

    % --- Plot 2: Beta' convergence with sum ---
    fig2 = figure('Position', [100, 50, fig_width, fig_height]);
    plot(1:k, history.beta_prime', 'LineWidth', 1.5);
    hold on;
    plot(1:k, beta_prime_sum_history, 'k-', 'LineWidth', 2);
    yline(beta_prime_sym, 'r--', 'LineWidth', 1.5);
    yline(beta_prime_sum_opt, 'k--', 'LineWidth', 1.5);
    xlabel('Iteration k');
    ylabel('\beta_i''(k)');
    title('Transformed Risk Parameter Convergence');
    legend_str = cell(n+3, 1);
    for i = 1:n
        legend_str{i} = sprintf('Operator %d', i);
    end
    legend_str{n+1} = '\Sigma\beta_i''';
    legend_str{n+2} = sprintf('\\beta_{sym}^* (%.4f)', beta_prime_sym);
    legend_str{n+3} = sprintf('\\Sigma\\beta^* (%.4f)', beta_prime_sum_opt);
    legend(legend_str, 'Location', 'best');
    grid on;
    if save_figs
        saveas(fig2, fullfile('figs', ['beta_convergence_' config_str '.png']));
    end

    % --- Plot 3: f and S convergence ---
    fig3 = figure('Position', [150, 50, fig_width, fig_height]);
    yyaxis left;
    plot(1:k, history.f, 'b-', 'LineWidth', 1.5);
    hold on;
    yline(0.5, 'b--', 'LineWidth', 1);
    ylabel('f(k)');
    yyaxis right;
    semilogy(1:k, abs(history.S), 'r-', 'LineWidth', 1.5);
    hold on;
    semilogy(1:k, history.alpha_change, 'm-', 'LineWidth', 1.5);
    yline(tau, 'r--', 'LineWidth', 1);
    yline(alpha_tol, 'm--', 'LineWidth', 1);
    ylabel('|S(k)|, |\Delta\alpha|');
    xlabel('Iteration k');
    title('Convergence Metrics');
    legend('f', 'f = 0.5', '|S|', '|\Delta\alpha|', '\tau', '\alpha_{tol}', 'Location', 'best');
    grid on;
    if save_figs
        saveas(fig3, fullfile('figs', ['convergence_metrics_' config_str '.png']));
    end

    % --- Plot 4: Individual profits ---
    fig4 = figure('Position', [200, 50, fig_width, fig_height]);
    plot(1:k, history.Pi', 'LineWidth', 1.5);
    xlabel('Iteration k');
    ylabel('\Pi_i(k)');
    title('Individual Profit Convergence');
    legend_str = cell(n, 1);
    for i = 1:n
        legend_str{i} = sprintf('Operator %d', i);
    end
    legend(legend_str, 'Location', 'best');
    grid on;
    if save_figs
        saveas(fig4, fullfile('figs', ['individual_profit_' config_str '.png']));
    end

    % --- Plot 5: Individual objectives ---
    fig5 = figure('Position', [250, 50, fig_width, fig_height]);
    plot(1:k, history.J', 'LineWidth', 1.5);
    xlabel('Iteration k');
    ylabel('J_i(k)');
    title('Individual Objective Convergence');
    legend(legend_str, 'Location', 'best');
    grid on;
    if save_figs
        saveas(fig5, fullfile('figs', ['individual_objective_' config_str '.png']));
    end

    % --- Plot 6: Aggregate metrics ---
    fig6 = figure('Position', [300, 50, fig_width, fig_height]);
    g_history = history.f .* (1 - history.f) * xr0^2;
    plot(1:k, sum(history.Pi, 1), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(1:k, g_history, 'r-', 'LineWidth', 1.5);
    yline(0.25*xr0^2, 'r--', 'LineWidth', 1);
    yline(sum(Pi_sym), 'b--', 'LineWidth', 1);
    xlabel('Iteration k');
    ylabel('Value');
    title('Aggregate Metrics Convergence');
    legend('\Sigma\Pi_i', 'g = f(1-f)x_r^{o2}', 'g_{max}', '\Pi_{sym}', ...
        'Location', 'best');
    grid on;
    if save_figs
        saveas(fig6, fullfile('figs', ['aggregate_metrics_' config_str '.png']));
    end

    if save_figs
        fprintf('Figures saved to figs/ folder with prefix: %s\n', config_str);
    end
end

%% ========================================================================
%                         HELPER FUNCTIONS
% =========================================================================

function z_bar = run_consensus(z0, L, epsilon, max_iter, tol)
    % Run consensus protocol until convergence
    % Returns the average (consensus value)
    z = z0;
    for iter = 1:max_iter
        z_new = z - epsilon * L * z;
        if norm(z_new - z) < tol
            break;
        end
        z = z_new;
    end
    z_bar = mean(z);
end

function [alpha_star, f, Pi, J] = compute_nash_centralized(beta_prime, ...
    n, a2, x0, xL0, xr0, eta, gamma, sigma_r2, sigma_L2, sigma_rL)
    % Compute Nash equilibrium for given beta' using Theorem 1
    % Inputs:
    %   beta_prime: n x 1 vector of transformed risk parameters
    %   Other parameters as defined in main script
    % Outputs:
    %   alpha_star: Nash equilibrium bidding fractions
    %   f: aggregate bidding fraction
    %   Pi: individual profits
    %   J: individual objectives

    % Build matrix A and vector B from Theorem 1 (Eq. 3.5)
    A = zeros(n, n);
    B = zeros(n, 1);

    for i = 1:n
        bp_i = beta_prime(i);
        for j = 1:n
            if i == j
                A(i, j) = 2 * (1 + 2*bp_i);
            else
                A(i, j) = (1 + 4*bp_i) * gamma(i, j);
            end
        end
        B(i) = (1 + 4*bp_i) * sum(gamma(i, :));
    end

    % Solve for Nash equilibrium
    alpha_star = A \ B;

    % Clamp to valid range (should be in [0,1] for valid beta')
    alpha_star = max(0, min(1, alpha_star));

    % Compute f
    f = sum(alpha_star .* eta);

    % Compute profits using Eq. (3.17)
    Pi = compute_profit(alpha_star, n, 1, 1, x0, xL0, xr0, eta, ...
        sigma_r2, sigma_L2, sigma_rL);
    % Scale by a2 (the full formula has 2*a2 coefficient)
    Pi = 2*a2 * alpha_star .* x0 .* xr0 .* (1 - f) + x0 + ...
         2*a2 * x0 .* (xL0 - xr0);

    % Add uncertainty terms if nonzero
    if sigma_r2 > 0 || sigma_L2 > 0 || sigma_rL ~= 0
        % Simplified: add covariance contribution
        Pi = Pi + 2*a2 * (sigma_rL - sigma_r2/n) * ones(n,1);
    end

    % Compute objectives J_i = Pi_i - beta_i' * 4*a2 * (x_r^o)^2 * (1-f)^2
    % (Simplified formula for E[(lambda - lambda^o)^2])
    price_var = 4 * a2^2 * (xr0^2 * (1-f)^2 + sigma_L2 + sigma_r2 - 2*sigma_rL);
    J = Pi - (beta_prime / a2) .* price_var;
end

function Pi = compute_profit(alpha, n, a2, a1, x0, xL0, xr0, eta, ...
    sigma_r2, sigma_L2, sigma_rL)
    % Compute individual profits given bidding fractions
    % Based on Eq. (3.17)

    f = sum(alpha .* eta);

    % Main profit term: 2*a2 * alpha_i * x_i^o * x_r^o * (1-f)
    Pi = 2*a2 * alpha .* x0 .* xr0 .* (1 - f);

    % Add linear term: a1 * x_i^o
    Pi = Pi + a1 * x0;

    % Add net load term: 2*a2 * x_i^o * (x_L^o - x_r^o)
    Pi = Pi + 2*a2 * x0 .* (xL0 - xr0);

    % Uncertainty contribution (if any)
    if sigma_rL ~= 0 || sigma_r2 > 0
        % Simplified model: shared equally among operators
        Pi = Pi + 2*a2 * (sigma_rL - sigma_r2/n) * ones(n, 1);
    end
end
