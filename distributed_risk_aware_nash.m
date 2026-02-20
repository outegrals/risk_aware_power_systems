%% Distributed Risk-Aware Nash Optimization Algorithm (IMPROVED VERSION)
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
% SINGLE TEST CASE VERSION - for debugging and parameter exploration
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

% Initial transformed risk parameters beta_i' = a2 * beta_i (user choice)
% Negative = risk-seeking, Zero = risk-neutral, Positive = risk-averse
% TEST CASE V (Risk-averse) - from paper
beta_prime_init = [-.23; -.1; -.21];

% Communication graph (adjacency matrix) - must be connected
% Ring topology: 1 -- 2 -- 3 -- 1
Adj = [0 1 1;
       1 0 1;
       1 1 0];  % Fully connected

% Algorithm parameters
epsilon_consensus = 0.2;  % Consensus step size (< 1/max_degree)
epsilon2 = 0.5;          % Beta update step size (NEW: increased for faster convergence)
tau = 1e-6;              % Stopping tolerance for |S| = |f - 0.5|
beta_tol = 1e-6;         % Stopping tolerance for beta convergence
min_iter = 10;           % Minimum iterations before checking convergence
max_outer_iter = 1000;   % Max iterations for outer loop
max_consensus_iter = 100;  % Max iterations for consensus
consensus_tol = 1e-10;   % Consensus convergence tolerance
delta = 0.01;            % Safety margin for beta' > -0.25

% Visualization parameters
save_figs = false;       % Set to true to save figures to figs/ folder
use_subplots = true;     % Set to true for one big subplot, false for individual figures

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
fprintf('   IMPROVED DISTRIBUTED RISK-AWARE NASH OPTIMIZATION\n');
fprintf('   (Direct Alpha Computation - No Damped Iteration)\n');
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
%                    PART 2: IMPROVED DISTRIBUTED ALGORITHM
% =========================================================================

fprintf('=========================================================\n');
fprintf('         PART 2: IMPROVED DISTRIBUTED ALGORITHM\n');
fprintf('      (Direct Alpha Computation - No Iteration)\n');
fprintf('=========================================================\n\n');

% Initialize variables
beta_prime = beta_prime_init;

% History for plotting
history.alpha = zeros(n, max_outer_iter);
history.beta_prime = zeros(n, max_outer_iter);
history.f = zeros(1, max_outer_iter);
history.S = zeros(1, max_outer_iter);
history.Pi = zeros(n, max_outer_iter);
history.J = zeros(n, max_outer_iter);
history.sum_beta_prime = zeros(1, max_outer_iter);

% --- Step 1: One-time consensus for aggregate quantities ---
fprintf('Step 1: One-time consensus for n, x_r^o, eta_i\n');

% Compute n via consensus using z_n
z_n = zeros(n, 1); z_n(1) = 1;
z_n_bar = run_consensus(z_n, L, epsilon_consensus, max_consensus_iter, consensus_tol);
n_est = 1 / z_n_bar;
fprintf('  Estimated n = %.4f (true = %d)\n', n_est, n);

% Compute x_r^o via consensus using z_xr0
z_xr0 = x0;
z_xr0_bar = run_consensus(z_xr0, L, epsilon_consensus, max_consensus_iter, consensus_tol);
xr0_est = n_est * z_xr0_bar;
fprintf('  Estimated x_r^o = %.4f (true = %.2f)\n', xr0_est, xr0);

% Compute eta_i locally using z_eta
z_eta = x0 / xr0_est;
eta_est = z_eta;
fprintf('  Estimated eta = [%s]\n\n', num2str(eta_est', '%.4f '));

% --- Step 2: Iterative optimization (IMPROVED) ---
fprintf('Step 2: Iterative optimization loop (IMPROVED ALGORITHM)\n');
fprintf('  Initial beta'' = [%s]\n\n', num2str(beta_prime', '%.4f '));

converged = false;
sum_beta_prev = sum(beta_prime);

fprintf('Iter |  sum(beta'')  |    f(alpha)   |      |S|      |  Converged?\n');
fprintf('-----+---------------+---------------+---------------+-------------\n');

for k = 1:max_outer_iter
    % Store history
    history.beta_prime(:, k) = beta_prime;
    
    % ====================================================================
    % STEP 1: CONSENSUS ON SUM OF BETA' (NEW: only one consensus round!)
    % Using z_beta(k) as consensus variable at iteration k
    % ====================================================================
    z_beta = beta_prime;
    z_beta_bar = run_consensus(z_beta, L, epsilon_consensus, max_consensus_iter, consensus_tol);
    z_beta_sum = n_est * z_beta_bar;  % Consensus estimate of sum(beta')
    sum_beta_prime = z_beta_sum;
    history.sum_beta_prime(k) = sum_beta_prime;

    % ====================================================================
    % STEP 2: DIRECT COMPUTATION OF NASH EQUILIBRIUM (NEW: no iteration!)
    % From Lemma 2: alpha*_i = (1+4*beta'_i)/eta_i * [1/(1+n+4*sum(beta'_j))]
    % Using consensus terms: z_eta and z_beta_sum(k)
    % ====================================================================
    denominator = 1 + n_est + 4 * z_beta_sum;
    alpha = (1 + 4 * beta_prime) ./ z_eta / denominator;

    % Clamp to [0,1] for safety (should be satisfied if beta' > -0.25)
    alpha = max(0, min(1, alpha));
    history.alpha(:, k) = alpha;

    % ====================================================================
    % STEP 3: COMPUTE f FROM ACTUAL (CLAMPED) ALPHA
    % NOTE: Must use actual alpha, not formula, because clamping
    % breaks the theoretical relationship f = (n+4*sum(beta'))/(1+n+4*sum(beta'))
    % Using consensus term z_eta
    % ====================================================================
    f = sum(alpha .* z_eta);  % Actual f from clamped alpha using z_eta
    history.f(k) = f;
    
    % ====================================================================
    % STEP 4: COMPUTE ERROR SIGNAL
    % ====================================================================
    S = f - 0.5;
    history.S(k) = S;
    
    % ====================================================================
    % STEP 5: UPDATE BETA' (gradient descent on Pareto error)
    % FIX: Iterative deficit redistribution to ensure total change
    % in sum(beta') matches the intended n*epsilon2*S
    % ====================================================================
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

    % Project to feasible region: beta' > -0.25 + delta
    beta_prime = max(beta_prime, lower_bound);
    
    % Compute profits and objectives at current state
    [~, ~, Pi_k, J_k] = compute_nash_centralized(beta_prime, n, a2, ...
        x0, xL0, xr0, eta, gamma, sigma_r2, sigma_L2, sigma_rL);
    Pi_actual = compute_profit(alpha, n, a2, a1, x0, xL0, xr0, eta, ...
        sigma_r2, sigma_L2, sigma_rL);
    price_var = 4 * a2^2 * (xr0^2 * (1-f)^2 + sigma_L2 + sigma_r2 - 2*sigma_rL);
    J_actual = Pi_actual - (beta_prime / a2) .* price_var;
    
    history.Pi(:, k) = Pi_actual;
    history.J(:, k) = J_actual;
    
    % Print progress every 10 iterations or if converged
    if mod(k, 10) == 0 || (abs(S) < tau && abs(sum_beta_prime - sum_beta_prev) < beta_tol && k > min_iter)
        fprintf('%4d | %12.6f | %12.6f | %13.6e |', ...
            k, sum_beta_prime, f, abs(S));
    end
    
    % Check convergence
    if abs(S) < tau && abs(sum_beta_prime - sum_beta_prev) < beta_tol && k > min_iter
        converged = true;
        fprintf('     YES\n');
        fprintf('\n*** CONVERGED at iteration %d ***\n\n', k);
        break;
    else
        if mod(k, 10) == 0
            fprintf('      NO\n');
        end
    end
    
    sum_beta_prev = sum_beta_prime;
end

if ~converged
    fprintf('\n*** WARNING: Did not converge within %d iterations ***\n\n', max_outer_iter);
end

% Trim history
history.alpha = history.alpha(:, 1:k);
history.beta_prime = history.beta_prime(:, 1:k);
history.f = history.f(1:k);
history.S = history.S(1:k);
history.Pi = history.Pi(:, 1:k);
history.J = history.J(:, 1:k);
history.sum_beta_prime = history.sum_beta_prime(1:k);

%% ========================================================================
%                         FINAL RESULTS
% =========================================================================

fprintf('=========================================================\n');
fprintf('                    FINAL RESULTS\n');
fprintf('=========================================================\n\n');

fprintf('Algorithm: IMPROVED (Direct Alpha Computation)\n');
fprintf('Iterations: %d\n', k);
fprintf('Convergence: %s\n\n', mat2str(converged));

fprintf('Final bidding fractions alpha:\n');
fprintf('  [%s]\n\n', num2str(alpha', '%.6f '));

fprintf('Final transformed risk parameters beta'':\n');
fprintf('  [%s]\n', num2str(beta_prime', '%.6f '));
fprintf('  sum(beta'') = %.6f (target: %.6f)\n\n', sum(beta_prime), beta_prime_sum_target);

fprintf('Final aggregate bidding fraction f:\n');
fprintf('  f = %.6f (target: 0.5)\n', f);
fprintf('  |S| = %.6e\n\n', abs(S));

fprintf('Market efficiency g = f*(1-f)*x_r^o^2:\n');
fprintf('  g = %.4f\n', f*(1-f)*xr0^2);
fprintf('  g_max = %.4f (Pareto optimal)\n', 0.25*xr0^2);
fprintf('  Efficiency: %.2f%%\n\n', 100*f*(1-f)*xr0^2/(0.25*xr0^2));

fprintf('Individual profits Pi_i:\n');
for i = 1:n
    fprintf('  Operator %d: %.4f\n', i, history.Pi(i, end));
end
fprintf('  Aggregate: %.4f\n', sum(history.Pi(:, end)));
fprintf('  Target (Pareto): %.4f\n\n', sum(Pi_sym));

fprintf('Comparison with Risk-Neutral Nash:\n');
fprintf('  Profit improvement: %.2f%%\n', ...
    100*(sum(history.Pi(:,end)) - sum(Pi_RN))/sum(Pi_RN));

%% ========================================================================
%                         VISUALIZATION
% =========================================================================

beta_prime_sum_history = sum(history.beta_prime, 1);
beta_prime_sum_opt = beta_prime_sum_target;

fig_width = 800;
fig_height = 500;

if use_subplots
    % === SUBPLOT MODE ===
    figure('Position', [50, 50, 1400, 900]);
    sgtitle(sprintf('Improved Distributed Risk-Aware Nash Optimization (Test Case: beta''_0=[%s])', ...
        num2str(beta_prime_init', '%.1f ')), 'FontSize', 14);
    
    % Plot 1: Sum of beta'
    subplot(2,3,1);
    plot(1:k, beta_prime_sum_history, 'b-', 'LineWidth', 2);
    hold on;
    yline(beta_prime_sum_opt, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Target');
    xlabel('Iteration k');
    ylabel('\Sigma\beta''_i(k)');
    title('Sum of Risk Parameters');
    legend('Location', 'best');
    grid on;
    
    % Plot 2: Individual beta'
    subplot(2,3,2);
    plot(1:k, history.beta_prime', 'LineWidth', 1.5);
    hold on;
    yline(beta_prime_sym, 'k--', 'LineWidth', 1);
    xlabel('Iteration k');
    ylabel('\beta''_i(k)');
    title('Individual Risk Parameters');
    legend_str = cell(n+1, 1);
    for i = 1:n
        legend_str{i} = sprintf('Operator %d', i);
    end
    legend_str{n+1} = sprintf('\\beta_{sym} = %.4f', beta_prime_sym);
    legend(legend_str, 'Location', 'best');
    grid on;
    
    % Plot 3: f convergence
    subplot(2,3,3);
    plot(1:k, history.f, 'b-', 'LineWidth', 2);
    hold on;
    yline(0.5, 'k--', 'LineWidth', 1.5);
    xlabel('Iteration k');
    ylabel('f(k)');
    title('Aggregate Bidding Fraction');
    legend('f(k)', 'Target = 0.5', 'Location', 'best');
    grid on;
    
    % Plot 4: Individual alpha
    subplot(2,3,4);
    plot(1:k, history.alpha', 'LineWidth', 1.5);
    xlabel('Iteration k');
    ylabel('\alpha_i(k)');
    title('Individual Bidding Fractions (Direct Computation - No Oscillations!)');
    legend_str = cell(n, 1);
    for i = 1:n
        legend_str{i} = sprintf('Operator %d', i);
    end
    legend(legend_str, 'Location', 'best');
    grid on;
    
    % Plot 5: Error signal |S|
    subplot(2,3,5);
    semilogy(1:k, abs(history.S), 'r-', 'LineWidth', 2);
    hold on;
    yline(tau, 'k--', 'LineWidth', 1);
    xlabel('Iteration k');
    ylabel('|S(k)| = |f - 0.5|');
    title('Pareto Error Signal (log scale)');
    legend('|S(k)|', '\tau', 'Location', 'best');
    grid on;
    
    % Plot 6: Aggregate profit and efficiency
    subplot(2,3,6);
    g_history = history.f .* (1 - history.f) * xr0^2;
    plot(1:k, sum(history.Pi, 1), 'b-', 'LineWidth', 2);
    hold on;
    plot(1:k, g_history, 'r-', 'LineWidth', 2);
    yline(0.25*xr0^2, 'r--', 'LineWidth', 1);
    yline(sum(Pi_sym), 'b--', 'LineWidth', 1);
    xlabel('Iteration k');
    ylabel('Value');
    title('Aggregate Metrics');
    legend('\Sigma\Pi_i', 'g = f(1-f)x_r^{o2}', 'g_{max}', '\Pi_{Pareto}', ...
        'Location', 'best');
    grid on;
    
    if save_figs
        saveas(gcf, 'figs/improved_algorithm_results.png');
        fprintf('Figure saved to: figs/improved_algorithm_results.png\n');
    end
    
else
    % === INDIVIDUAL FIGURE MODE ===
    % (Similar to original, omitted for brevity)
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