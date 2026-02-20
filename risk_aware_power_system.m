%% Risk-Aware Power System - Multiple Test Cases Comparison (IMPROVED)
% Implementation of IMPROVED Algorithm from the paper
% "Distributed Risk-Aware Bidding Strategy for Incorporating Renewable
% Generation into Real-Time Electricity Market"
%
% KEY IMPROVEMENT: Direct computation of alpha* instead of damped iteration
%   - Eliminates oscillations
%   - Faster convergence
%   - Only one consensus round per iteration (instead of two)
%   - Single tuning parameter (epsilon2) instead of two
%
% This version runs MULTIPLE TEST CASES and displays results on combined plots.
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

% Communication graph (adjacency matrix) - must be connected
Adj = [0 1 1;
       1 0 1;
       1 1 0];  % Fully connected

% Algorithm parameters
epsilon_consensus = 0.2;  % Consensus step size (< 1/max_degree)
epsilon2 = 0.5;          % Beta update step size
tau = 1e-6;              % Stopping tolerance for |S| = |f - 0.5|
beta_tol = 1e-6;         % Stopping tolerance for beta convergence
min_iter = 10;           % Minimum iterations before checking convergence
max_outer_iter = 1000;   % Max iterations for outer loop
max_consensus_iter = 100;  % Max iterations for consensus
consensus_tol = 1e-10;   % Consensus convergence tolerance
delta = 0.01;            % Safety margin for beta' > -0.25

% Visualization parameters
save_figs = true;           % Set to true to save figures to figs/ folder

%% ========================================================================
%                         TEST CASES DEFINITION
% =========================================================================
% Each row is a test case: [beta1', beta2', beta3']
% Negative = risk-seeking, Zero = risk-neutral, Positive = risk-averse

test_cases = {
    'Case I',   [0; 0; 0],           'Risk-neutral';
    'Case II',  [-1; -1; -1],        'Risk-seeking';
    'Case III', [1; 1; 1],           'Risk-averse';
    'Case IV',  [-1; 1; 1],          'Mixed';
    'Case V',   [-0.25; 0; 0.2],    'Mixed (mild)';
    'Case VI',   [-0.25; 0; 0],    'Mixed (mild)';
};

num_cases = size(test_cases, 1);
case_names = test_cases(:, 1);
case_betas = test_cases(:, 2);
case_profiles = test_cases(:, 3);

%% ========================================================================
%                         DERIVED QUANTITIES
% =========================================================================

xr0 = sum(x0);           % Aggregate renewable forecast
eta = x0 / xr0;          % Forecast ratios eta_i = x_i^o / x_r^o
gamma = x0' ./ x0;       % gamma_{ij} = x_j^o / x_i^o (n x n matrix)

% Laplacian matrix for consensus
Deg = diag(sum(Adj, 2));
L = Deg - Adj;

% Storage for results from all test cases
all_results = cell(num_cases, 1);

fprintf('=========================================================\n');
fprintf('   DISTRIBUTED RISK-AWARE NASH OPTIMIZATION\n');
fprintf('         MULTIPLE TEST CASES COMPARISON\n');
fprintf('=========================================================\n\n');
fprintf('Problem Parameters:\n');
fprintf('  n = %d operators\n', n);
fprintf('  a2 = %.2f, a1 = %.2f, a0 = %.2f\n', a2, a1, a0);
fprintf('  x0 = [%s] MW\n', num2str(x0'));
fprintf('  xL0 = %.2f MW\n', xL0);
fprintf('  x_r^o = %.2f MW\n', xr0);
fprintf('  eta = [%s]\n\n', num2str(eta', '%.4f '));

fprintf('Test Cases:\n');
for tc = 1:num_cases
    fprintf('  %s: beta'' = [%s] (%s)\n', case_names{tc}, ...
        num2str(case_betas{tc}', '%.2f '), case_profiles{tc});
end
fprintf('\n');

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
beta_prime_sum_opt = n * beta_prime_sym;  % Optimal sum of beta'

fprintf('Theoretical beta_sym'' = (1-n)/(4n) = %.6f\n', beta_prime_sym);
fprintf('Theoretical sum(beta'') = (1-n)/4 = %.6f\n', beta_prime_sum_opt);
fprintf('Theoretical beta_sym = beta_sym''/a2 = %.6f\n', beta_sym);

[alpha_sym, f_sym, Pi_sym, J_sym] = compute_nash_centralized(...
    beta_prime_sym * ones(n,1), n, a2, x0, xL0, xr0, eta, gamma, ...
    sigma_r2, sigma_L2, sigma_rL);

fprintf('alpha_sym = [%s]\n', num2str(alpha_sym', '%.6f '));
fprintf('f_sym = %.6f (target: 0.5)\n', f_sym);
fprintf('Aggregate profit Pi_sym = %.4f\n', sum(Pi_sym));
fprintf('g_sym = f*(1-f)*x_r^o^2 = %.4f\n\n', f_sym*(1-f_sym)*xr0^2);

%% ========================================================================
%                    PART 2: DISTRIBUTED ALGORITHM (ALL TEST CASES)
% =========================================================================

fprintf('=========================================================\n');
fprintf('      PART 2: IMPROVED DISTRIBUTED ALGORITHM\n');
fprintf('      (Direct Alpha Computation - No Iteration)\n');
fprintf('=========================================================\n\n');

% --- One-time consensus for aggregate quantities (same for all cases) ---
fprintf('One-time consensus for n, x_r^o, eta_i\n');

% Compute n via consensus (Eq. 5.3) using z_n
z_n = zeros(n, 1); z_n(1) = 1;
z_n_bar = run_consensus(z_n, L, epsilon_consensus, max_consensus_iter, consensus_tol);
n_est = 1 / z_n_bar;
fprintf('  Estimated n = %.4f (true = %d)\n', n_est, n);

% Compute x_r^o via consensus (Eq. 5.4) using z_xr0
z_xr0 = x0;
z_xr0_bar = run_consensus(z_xr0, L, epsilon_consensus, max_consensus_iter, consensus_tol);
xr0_est = n_est * z_xr0_bar;
fprintf('  Estimated x_r^o = %.4f (true = %.2f)\n', xr0_est, xr0);

% Compute eta_i locally (Eq. 5.5) using z_eta
z_eta = x0 / xr0_est;
eta_est = z_eta;
fprintf('  Estimated eta = [%s]\n\n', num2str(eta_est', '%.4f '));

% --- Run distributed algorithm for each test case ---
for tc = 1:num_cases
    fprintf('--- %s: %s ---\n', case_names{tc}, case_profiles{tc});

    % Initialize variables for this test case
    beta_prime_init = case_betas{tc};
    beta_prime = beta_prime_init;

    % Clamp initial beta' to feasible region (beta' > -0.25)
    beta_prime = max(beta_prime, -0.25 + delta);

    fprintf('  Initial beta'' = [%s]', num2str(beta_prime_init', '%.4f '));
    if any(beta_prime_init < -0.25 + delta)
        fprintf(' -> clamped to [%s]', num2str(beta_prime', '%.4f '));
    end
    fprintf('\n');

    % History for plotting
    history.alpha = zeros(n, max_outer_iter);
    history.beta_prime = zeros(n, max_outer_iter);
    history.f = zeros(1, max_outer_iter);
    history.S = zeros(1, max_outer_iter);
    history.Pi = zeros(n, max_outer_iter);
    history.J = zeros(n, max_outer_iter);
    history.sum_beta_prime = zeros(1, max_outer_iter);

    converged = false;
    sum_beta_prev = sum(beta_prime);

    for k = 1:max_outer_iter
        % Store history
        history.beta_prime(:, k) = beta_prime;

        % ============================================================
        % STEP 1: CONSENSUS ON SUM OF BETA' (one consensus round)
        % Using z_beta(k) as consensus variable at iteration k
        % ============================================================
        z_beta = beta_prime;
        z_beta_bar = run_consensus(z_beta, L, epsilon_consensus, max_consensus_iter, consensus_tol);
        z_beta_sum = n_est * z_beta_bar;  % Consensus estimate of sum(beta')
        sum_beta_prime = z_beta_sum;
        history.sum_beta_prime(k) = sum_beta_prime;

        % ============================================================
        % STEP 2: DIRECT COMPUTATION OF NASH EQUILIBRIUM
        % alpha*_i = (1+4*beta'_i)/eta_i * [1/(1+n+4*sum(beta'_j))]
        % Using consensus terms: z_eta and z_beta_sum(k)
        % ============================================================
        denominator = 1 + n_est + 4 * z_beta_sum;
        alpha = (1 + 4 * beta_prime) ./ z_eta / denominator;
        alpha = max(0, min(1, alpha));  % Clamp to [0,1]
        history.alpha(:, k) = alpha;

        % ============================================================
        % STEP 3: COMPUTE f FROM ACTUAL (CLAMPED) ALPHA
        % NOTE: Must use actual alpha, not formula, because clamping
        % breaks the theoretical relationship f = (n+4*sum(beta'))/(1+n+4*sum(beta'))
        % Using consensus term z_eta
        % ============================================================
        f = sum(alpha .* z_eta);  % Actual f from clamped alpha using z_eta
        history.f(k) = f;

        % ============================================================
        % STEP 4: COMPUTE ERROR SIGNAL
        % ============================================================
        S = f - 0.5;
        history.S(k) = S;

        % ============================================================
        % STEP 5: UPDATE BETA' (gradient descent on Pareto error)
        % FIX: Iterative deficit redistribution to ensure total change
        % in sum(beta') matches the intended n*epsilon2*S
        % ============================================================
        lower_bound = -0.25 + delta;

        if S > 0  % Need to decrease beta' values
            % Total desired decrease in sum(beta')
            remaining_delta = n * epsilon2 * S;

            % Iteratively distribute the change, accounting for clamping
            for redistrib_iter = 1:10
                % Find operators with headroom (not at bound)
                headroom = beta_prime - lower_bound;
                has_headroom = headroom > 1e-10;
                n_free = sum(has_headroom);

                if n_free == 0 || remaining_delta < 1e-12
                    break;
                end

                % Distribute remaining delta among free operators
                per_operator = remaining_delta / n_free;

                % Compute actual change (limited by headroom)
                actual_change = min(per_operator, headroom) .* has_headroom;

                % Apply change
                beta_prime = beta_prime - actual_change;

                % Update remaining delta for next iteration
                remaining_delta = remaining_delta - sum(actual_change);
            end
        else  % S <= 0, need to increase beta' values (no upper bound)
            beta_prime = beta_prime - epsilon2 * S;
        end

        % Final safety clamping
        beta_prime = max(beta_prime, lower_bound);

        % Compute profits and objectives at current state
        Pi_actual = compute_profit(alpha, n, a2, a1, x0, xL0, xr0, eta, ...
            sigma_r2, sigma_L2, sigma_rL);
        price_var = 4 * a2^2 * (xr0^2 * (1-f)^2 + sigma_L2 + sigma_r2 - 2*sigma_rL);
        J_actual = Pi_actual - (beta_prime / a2) .* price_var;
        history.Pi(:, k) = Pi_actual;
        history.J(:, k) = J_actual;

        % Check convergence
        if k >= min_iter && abs(S) < tau && abs(sum_beta_prime - sum_beta_prev) < beta_tol
            fprintf('  Converged at iteration %d: |S| = %.2e\n', k, abs(S));
            converged = true;
            break;
        end

        sum_beta_prev = sum_beta_prime;
    end

    if ~converged
        fprintf('  Did not converge within %d iterations\n', max_outer_iter);
    end

    % Trim history and store results
    history.alpha = history.alpha(:, 1:k);
    history.beta_prime = history.beta_prime(:, 1:k);
    history.f = history.f(1:k);
    history.S = history.S(1:k);
    history.Pi = history.Pi(:, 1:k);
    history.J = history.J(:, 1:k);
    history.sum_beta_prime = history.sum_beta_prime(1:k);
    history.iterations = k;
    history.converged = converged;
    history.beta_prime_init = beta_prime_init;

    % Store in results array
    all_results{tc} = history;

    fprintf('  Final f = %.6f, sum(Pi) = %.4f\n\n', history.f(end), sum(history.Pi(:, end)));
end

%% ========================================================================
%                         FINAL RESULTS SUMMARY
% =========================================================================

fprintf('\n=========================================================\n');
fprintf('                    FINAL RESULTS SUMMARY\n');
fprintf('=========================================================\n\n');

fprintf('Centralized Reference Values:\n');
fprintf('  Symmetric beta''_sym = %.6f\n', beta_prime_sym);
fprintf('  Optimal sum(beta'') = %.6f\n', beta_prime_sum_opt);
fprintf('  Symmetric f = %.6f\n', f_sym);
fprintf('  Symmetric sum(Pi) = %.4f\n', sum(Pi_sym));
fprintf('  Maximum g = %.4f\n\n', 0.25*xr0^2);

fprintf('%-10s | %6s | %10s | %10s | %10s | %10s | %8s\n', ...
    'Case', 'Iters', 'Final f', 'sum(beta'')', 'dist to opt', 'sum(Pi)', 'g');
fprintf('%s\n', repmat('-', 1, 82));

for tc = 1:num_cases
    res = all_results{tc};
    f_final = res.f(end);
    beta_sum = sum(res.beta_prime(:, end));
    dist_to_opt = abs(beta_sum - beta_prime_sum_opt);
    Pi_sum = sum(res.Pi(:, end));
    g_final = f_final * (1 - f_final) * xr0^2;

    fprintf('%-10s | %6d | %10.6f | %10.6f | %10.6f | %10.4f | %8.4f\n', ...
        case_names{tc}, res.iterations, f_final, beta_sum, dist_to_opt, Pi_sum, g_final);
end
fprintf('\n');

%% ========================================================================
%                         VISUALIZATION (COMBINED PLOTS)
% =========================================================================

% Create figs folder if saving is enabled
if save_figs && ~exist('figs', 'dir')
    mkdir('figs');
end

% Find max iterations across all cases for consistent x-axis
max_iters = max(cellfun(@(r) r.iterations, all_results));

% Color scheme for test cases
colors = lines(num_cases);

% Line styles for operators within each case
line_styles = {'-', '--', ':'};

fprintf('=========================================================\n');
fprintf('                    VISUALIZATION\n');
fprintf('=========================================================\n\n');

% === Figure 1: Sum of Beta' Convergence (all cases) ===
fig1 = figure('Position', [50, 50, 800, 500], 'Name', 'Sum Beta Convergence');
hold on;
for tc = 1:num_cases
    res = all_results{tc};
    beta_sum = sum(res.beta_prime, 1);
    plot(1:res.iterations, beta_sum, 'Color', colors(tc,:), ...
        'LineWidth', 2, 'DisplayName', case_names{tc});
end
yline(beta_prime_sum_opt, 'k--', 'LineWidth', 2, 'DisplayName', '\Sigma\beta''_{opt}');
xlabel('Iteration k');
ylabel('\Sigma\beta_i''(k)');
title('Sum of Risk Parameters Convergence (All Cases)');
legend('Location', 'best');
grid on;
if save_figs
    saveas(fig1, fullfile('figs', 'sum_beta_convergence_all_cases.png'));
end

% === Figure 2: f Convergence (all cases) ===
fig2 = figure('Position', [100, 50, 800, 500], 'Name', 'f Convergence');
hold on;
for tc = 1:num_cases
    res = all_results{tc};
    plot(1:res.iterations, res.f, 'Color', colors(tc,:), ...
        'LineWidth', 2, 'DisplayName', case_names{tc});
end
yline(0.5, 'k--', 'LineWidth', 2, 'DisplayName', 'f_{opt} = 0.5');
xlabel('Iteration k');
ylabel('f(k)');
title('Aggregate Bidding Fraction Convergence (All Cases)');
legend('Location', 'best');
grid on;
if save_figs
    saveas(fig2, fullfile('figs', 'f_convergence_all_cases.png'));
end

% === Figure 3: |S| Convergence (all cases) ===
fig3 = figure('Position', [150, 50, 800, 500], 'Name', 'S Convergence');
hold on;
for tc = 1:num_cases
    res = all_results{tc};
    semilogy(1:res.iterations, abs(res.S), 'Color', colors(tc,:), ...
        'LineWidth', 2, 'DisplayName', case_names{tc});
end
yline(tau, 'k--', 'LineWidth', 2, 'DisplayName', '\tau');
xlabel('Iteration k');
ylabel('|S(k)|');
title('Error Signal |S| = |f - 0.5| Convergence (All Cases)');
legend('Location', 'best');
grid on;
if save_figs
    saveas(fig3, fullfile('figs', 'S_convergence_all_cases.png'));
end

% === Figure 4: Aggregate Profit Convergence (all cases) ===
fig4 = figure('Position', [200, 50, 800, 500], 'Name', 'Profit Convergence');
hold on;
for tc = 1:num_cases
    res = all_results{tc};
    Pi_sum = sum(res.Pi, 1);
    plot(1:res.iterations, Pi_sum, 'Color', colors(tc,:), ...
        'LineWidth', 2, 'DisplayName', case_names{tc});
end
yline(sum(Pi_sym), 'k--', 'LineWidth', 2, 'DisplayName', '\Pi_{sym}');
xlabel('Iteration k');
ylabel('\Sigma\Pi_i(k)');
title('Aggregate Profit Convergence (All Cases)');
legend('Location', 'best');
grid on;
if save_figs
    saveas(fig4, fullfile('figs', 'profit_convergence_all_cases.png'));
end

% === Figure 5: g(alpha) Convergence (all cases) ===
fig5 = figure('Position', [250, 50, 800, 500], 'Name', 'g Convergence');
hold on;
for tc = 1:num_cases
    res = all_results{tc};
    g_hist = res.f .* (1 - res.f) * xr0^2;
    plot(1:res.iterations, g_hist, 'Color', colors(tc,:), ...
        'LineWidth', 2, 'DisplayName', case_names{tc});
end
yline(0.25*xr0^2, 'k--', 'LineWidth', 2, 'DisplayName', 'g_{max}');
xlabel('Iteration k');
ylabel('g(\alpha) = f(1-f)x_r^{o2}');
title('Market Efficiency g(\alpha) Convergence (All Cases)');
legend('Location', 'best');
grid on;
if save_figs
    saveas(fig5, fullfile('figs', 'g_convergence_all_cases.png'));
end

% === Figure 6: Individual Beta' for each case (subplots) ===
fig6 = figure('Position', [300, 50, 1200, 800], 'Name', 'Individual Beta Convergence');
for tc = 1:num_cases
    subplot(2, 3, tc);
    res = all_results{tc};
    plot(1:res.iterations, res.beta_prime', 'LineWidth', 1.5);
    hold on;
    beta_sum = sum(res.beta_prime, 1);
    plot(1:res.iterations, beta_sum, 'k-', 'LineWidth', 2);
    yline(beta_prime_sym, 'r--', 'LineWidth', 1.5);
    xlabel('Iteration k');
    ylabel('\beta_i''(k)');
    title(sprintf('%s (%s)', case_names{tc}, case_profiles{tc}));
    if tc == 1
        legend_str = cell(n+2, 1);
        for i = 1:n
            legend_str{i} = sprintf('Op %d', i);
        end
        legend_str{n+1} = '\Sigma\beta_i''';
        legend_str{n+2} = '\beta_{sym}''';
        legend(legend_str, 'Location', 'best', 'FontSize', 7);
    end
    grid on;
end
sgtitle('Individual \beta_i'' Convergence by Test Case', 'FontSize', 14, 'FontWeight', 'bold');
if save_figs
    saveas(fig6, fullfile('figs', 'beta_individual_all_cases.png'));
end

% === Figure 7: Individual Alpha for each case (subplots) ===
fig7 = figure('Position', [350, 50, 1200, 800], 'Name', 'Individual Alpha Convergence');
for tc = 1:num_cases
    subplot(2, 3, tc);
    res = all_results{tc};
    plot(1:res.iterations, res.alpha', 'LineWidth', 1.5);
    hold on;
    yline(0.5, 'k--', 'LineWidth', 1.5);
    xlabel('Iteration k');
    ylabel('\alpha_i(k)');
    title(sprintf('%s (%s)', case_names{tc}, case_profiles{tc}));
    if tc == 1
        legend_str = cell(n+1, 1);
        for i = 1:n
            legend_str{i} = sprintf('Op %d', i);
        end
        legend_str{n+1} = 'Target';
        legend(legend_str, 'Location', 'best', 'FontSize', 7);
    end
    grid on;
end
sgtitle('Individual \alpha_i Convergence by Test Case', 'FontSize', 14, 'FontWeight', 'bold');
if save_figs
    saveas(fig7, fullfile('figs', 'alpha_individual_all_cases.png'));
end

% === Figure 8: Individual Profit for each case (subplots) ===
fig8 = figure('Position', [400, 50, 1200, 800], 'Name', 'Individual Profit Convergence');
for tc = 1:num_cases
    subplot(2, 3, tc);
    res = all_results{tc};
    plot(1:res.iterations, res.Pi', 'LineWidth', 1.5);
    xlabel('Iteration k');
    ylabel('\Pi_i(k)');
    title(sprintf('%s (%s)', case_names{tc}, case_profiles{tc}));
    if tc == 1
        legend_str = cell(n, 1);
        for i = 1:n
            legend_str{i} = sprintf('Op %d', i);
        end
        legend(legend_str, 'Location', 'best', 'FontSize', 7);
    end
    grid on;
end
sgtitle('Individual Profit \Pi_i Convergence by Test Case', 'FontSize', 14, 'FontWeight', 'bold');
if save_figs
    saveas(fig8, fullfile('figs', 'profit_individual_all_cases.png'));
end

if save_figs
    fprintf('All figures saved to figs/ folder\n');
end

fprintf('Visualization complete.\n');

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
