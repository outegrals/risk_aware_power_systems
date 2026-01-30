% =========================================================================
% DISTRIBUTED RISK-AWARE NASH OPTIMIZATION FOR RENEWABLE ENERGY MARKETS
% =========================================================================
% This script implements Algorithms 1 and 2 from Section V
% Demonstrates convergence to Pareto-optimal Nash equilibrium
% =========================================================================

clear all; close all; clc;

%% ========================================================================
% CONFIGURATION
% =========================================================================

% Choose which algorithm to run
% 'iterative' - Algorithm 1: General distributed optimization
% 'symmetric' - Algorithm 2: Simplified symmetric case
ALGORITHM = 'iterative'; 

% Problem parameters
n = 3;                      % Number of operators
a2 = .8;                 % Cost coefficient ($/MW^2) - TYPICAL: 1e-4 to 1e-3
                           % NOTE: Large values (e.g., a2=0.8) are non-physical 
                           % and may cause unexpected behavior
a1 = 10;                   % Cost coefficient ($/MW)
a0 = 100;                  % Cost coefficient ($)

% Forecast generation (can modify for asymmetric case)
SYMMETRIC = true;          % Set to false for asymmetric forecasts
if SYMMETRIC
    x_o = 100 * ones(n,1); % Symmetric: all operators have same forecast (MW)
else
    x_o = [100; 80; 120];  % Asymmetric forecasts (MW)
end

x_r_o = sum(x_o);          % Total forecast
eta = x_o / x_r_o;         % Forecast fractions

% Display theoretical optimal beta for symmetric case
if SYMMETRIC || (max(x_o) - min(x_o)) < 1e-6
    beta_theory = (1 - n) / (4 * a2 * n);
    fprintf('\n==========================================================\n');
    fprintf('Theoretical symmetric optimal beta: %.6f\n', beta_theory);
    fprintf('==========================================================\n\n');
    
    if beta_theory > 0
        warning('Theoretical beta is POSITIVE (risk-averse). This is unusual and may indicate non-standard parameters.');
        fprintf('  For Pareto optimality, we typically expect beta < 0 (risk-seeking).\n');
        fprintf('  Your a2 = %.4f may be too large for typical power systems (usually 1e-4 to 1e-3).\n\n', a2);
    end
end

% Communication network (adjacency matrix)
% For now, use complete graph (all-to-all communication)
% Can modify for sparse networks
A_comm = ones(n,n) - eye(n); % Adjacency matrix (1 if connected, 0 otherwise)

% Algorithm parameters
epsilon1 = 0.1;            % Consensus step size
% Auto-scale epsilon2 based on a2 (smaller a2 needs smaller epsilon2)
epsilon2 = 0.01 / (4*a2*n); % Risk parameter update step size (auto-scaled)
tau = 1e-4;                % Convergence tolerance
max_iter = 500;            % Maximum iterations (increased for stability)

% Initial conditions
beta_init = zeros(n,1);    % Initial risk parameters (risk-neutral)
alpha_init = 0.5*ones(n,1);% Initial bidding fractions (will be overridden in algorithm)

%% ========================================================================
% RUN SELECTED ALGORITHM
% =========================================================================

switch ALGORITHM
    case 'symmetric'
        fprintf('Running Algorithm 2: Symmetric Pareto-Optimal Strategy\n');
        fprintf('==========================================================\n');
        [alpha_star, beta_star, ~] = algorithm_symmetric(n, a2, epsilon1);
        
        % For comparison, also compute risk-neutral Nash
        [alpha_RN, ~] = nash_equilibrium(zeros(n,1), x_o, a2);
        
        % Store for plotting (single iteration)
        alpha_history = alpha_star;
        beta_history = beta_star;
        f_history = sum(alpha_star .* eta);
        g_history = compute_g(alpha_star, eta, x_r_o);
        
        fprintf('\nResults:\n');
        fprintf('  Optimal beta: %.6f\n', beta_star(1));
        fprintf('  Optimal alpha: %.6f\n', alpha_star(1));
        fprintf('  f (should be 0.5): %.6f\n', f_history);
        fprintf('  g (aggregate profit factor): %.6f\n', g_history);
        
    case 'iterative'
        fprintf('Running Algorithm 1: Distributed Risk-Aware Nash Optimization\n');
        fprintf('==============================================================\n');
        [alpha_star, beta_star, history] = algorithm_iterative(...
            beta_init, alpha_init, x_o, a2, epsilon1, epsilon2, tau, max_iter, A_comm);
        
        % Extract history for plotting
        alpha_history = history.alpha;
        beta_history = history.beta;
        f_history = history.f;
        g_history = history.g;
        S_history = history.S;
        
        fprintf('\nConvergence achieved at iteration %d\n', length(f_history));
        fprintf('Final values:\n');
        fprintf('  f: %.6f (target: 0.5000)\n', f_history(end));
        fprintf('  g: %.6f (max: %.6f)\n', g_history(end), 0.25*x_r_o^2);
        fprintf('  Beta values:\n');
        for i = 1:n
            fprintf('    beta_%d: %.6f\n', i, beta_star(i));
        end
        
    otherwise
        error('Invalid algorithm selection');
end

%% ========================================================================
% COMPUTE BASELINE COMPARISONS
% =========================================================================

% Risk-neutral Nash equilibrium
beta_RN = zeros(n,1);
[alpha_RN, ~] = nash_equilibrium(beta_RN, x_o, a2);
f_RN = sum(alpha_RN .* eta);
g_RN = compute_g(alpha_RN, eta, x_r_o);

% Cooperative (centralized) solution
alpha_coop = 0.5 * ones(n,1);
beta_coop = zeros(n,1);
f_coop = sum(alpha_coop .* eta);
g_coop = compute_g(alpha_coop, eta, x_r_o);

% Pareto-optimal (from algorithm)
f_pareto = f_history(end);
g_pareto = g_history(end);

% Price of Anarchy calculations
PoA_RN = g_coop / g_RN;
PoA_pareto = g_coop / g_pareto;

%% ========================================================================
% ANALYTICAL OPTIMAL BETA CALCULATION
% =========================================================================

fprintf('\n==========================================================\n');
fprintf('ANALYTICAL OPTIMAL BETA CALCULATION\n');
fprintf('==========================================================\n');

% For symmetric case: beta_sym = (1-n)/(4*a2*n)
beta_analytical_sym = (1 - n) / (4 * a2 * n);
fprintf('Symmetric formula: beta = (1-n)/(4*a2*n)\n');
fprintf('  = (1-%d)/(4*%.4f*%d)\n', n, a2, n);
fprintf('  = %.6f (RISK-%s)\n\n', beta_analytical_sym, ...
    ternary(beta_analytical_sym < 0, 'SEEKING', 'AVERSE'));

% Compute Nash equilibrium with analytical beta
[alpha_analytical, ~] = nash_equilibrium(beta_analytical_sym * ones(n,1), x_o, a2);
f_analytical = sum(alpha_analytical .* eta);
g_analytical = compute_g(alpha_analytical, eta, x_r_o);

fprintf('With analytical beta:\n');
fprintf('  f = %.6f (target: 0.5000)\n', f_analytical);
fprintf('  g = %.6f (max: %.6f)\n', g_analytical, g_coop);
fprintf('  alpha values: [');
for i = 1:n
    fprintf('%.4f', alpha_analytical(i));
    if i < n, fprintf(', '); end
end
fprintf(']\n\n');

% For asymmetric case (duopoly constraint)
if n == 2
    beta_sum_target = -1/(4*a2);
    fprintf('Duopoly constraint: beta_1 + beta_2 = -1/(4*a2)\n');
    fprintf('  = -1/(4*%.4f) = %.6f\n', a2, beta_sum_target);
    fprintf('  Current sum: beta_1 + beta_2 = %.6f + %.6f = %.6f\n', ...
        beta_star(1), beta_star(2), sum(beta_star));
    fprintf('  Difference from target: %.6f\n', sum(beta_star) - beta_sum_target);
end

fprintf('==========================================================\n');

fprintf('\n==========================================================\n');
fprintf('COMPARISON OF STRATEGIES\n');
fprintf('==========================================================\n');
fprintf('Risk-Neutral Nash:\n');
fprintf('  f = %.4f, g = %.4f, PoA = %.4f\n', f_RN, g_RN, PoA_RN);
fprintf('Cooperative (Centralized):\n');
fprintf('  f = %.4f, g = %.4f, PoA = %.4f\n', f_coop, g_coop, 1.0);
fprintf('Pareto-Optimal (Distributed):\n');
fprintf('  f = %.4f, g = %.4f, PoA = %.4f\n', f_pareto, g_pareto, PoA_pareto);
fprintf('==========================================================\n');
fprintf('Improvement over Risk-Neutral: %.2f%%\n', 100*(g_pareto - g_RN)/g_RN);
fprintf('Gap to Cooperative: %.2f%%\n', 100*(g_coop - g_pareto)/g_coop);

%% ========================================================================
% VISUALIZATION
% =========================================================================

create_plots(ALGORITHM, alpha_history, beta_history, f_history, g_history, ...
             alpha_RN, g_RN, g_coop, n, eta, x_r_o);

%% ========================================================================
% FUNCTION DEFINITIONS
% =========================================================================

% -------------------------------------------------------------------------
function [alpha_star, beta_star, history] = algorithm_symmetric(n, a2, epsilon1)
% Algorithm 2: Symmetric Pareto-Optimal Strategy
    
    % Step 1: Compute n via consensus (simulated - all operators know n already)
    % In practice, this would use the consensus protocol
    n_computed = n;
    
    % Step 2: Compute optimal parameters
    beta_star = ((1 - n) / (4 * a2 * n)) * ones(n, 1);
    alpha_star = 0.5 * ones(n, 1);
    
    % Return empty history for symmetric case (non-iterative)
    history = struct();
end

% -------------------------------------------------------------------------
function [alpha_star, beta_star, history] = algorithm_iterative(...
    beta_init, alpha_init, x_o, a2, epsilon1, epsilon2, tau, max_iter, A_comm)
% Algorithm 1: Distributed Risk-Aware Nash Optimization

    n = length(x_o);
    x_r_o = sum(x_o);
    eta = x_o / x_r_o;
    
    % Initialize
    beta = beta_init;
    % Start from Nash equilibrium for initial beta (not arbitrary 0.5)
    [alpha, ~] = nash_equilibrium(beta_init, x_o, a2);
    k = 1;
    
    % Storage for history
    alpha_history = zeros(n, max_iter);
    beta_history = zeros(n, max_iter);
    f_history = zeros(max_iter, 1);
    g_history = zeros(max_iter, 1);
    S_history = zeros(max_iter, 1);
    
    % Iterative optimization
    converged = false;
    
    fprintf('\n--- Starting Iterations ---\n');
    fprintf('Initial state:\n');
    fprintf('  beta: [');
    for i = 1:n
        fprintf('%.4f', beta(i));
        if i < n, fprintf(', '); end
    end
    fprintf(']\n');
    fprintf('  alpha: [');
    for i = 1:n
        fprintf('%.4f', alpha(i));
        if i < n, fprintf(', '); end
    end
    fprintf(']\n');
    fprintf('  f = %.6f, target = 0.5000\n', sum(alpha .* eta));
    fprintf('  epsilon2 = %.6e\n\n', epsilon2);
    
    while ~converged && k <= max_iter
        % Step 1: Compute f(alpha(k)) via consensus
        f_k = consensus_sum(alpha .* x_o, A_comm, epsilon1) / x_r_o;
        
        % Step 2: Update alpha via Nash best response
        alpha = update_alpha(alpha, beta, f_k, eta, x_o, a2);
        
        % Step 3: Compute f(alpha(k+1)) via consensus
        f_k1 = consensus_sum(alpha .* x_o, A_comm, epsilon1) / x_r_o;
        
        % Step 4: Compute error signal S(k)
        S_k = f_k1 - 0.5;
        
        % Step 5: Update beta via gradient descent
        beta = beta - epsilon2 * S_k;
        
        % Enforce stability constraint: beta > -1/(4*a2) for all i
        % This is the theoretical minimum from the Nash equilibrium existence condition
        beta_min = -1/(4*a2) + 1e-6; % Small margin for numerical stability
        beta = max(beta, beta_min);
        
        % Store history
        alpha_history(:, k) = alpha;
        beta_history(:, k) = beta;
        f_history(k) = f_k1;
        g_history(k) = compute_g(alpha, eta, x_r_o);
        S_history(k) = S_k;
        
        % Check convergence
        if abs(S_k) < tau
            converged = true;
            fprintf('  Iteration %d: f = %.6f, |S| = %.6e (converged)\n', k, f_k1, abs(S_k));
            fprintf('    Final beta: [');
            for i = 1:n
                fprintf('%.4f', beta(i));
                if i < n, fprintf(', '); end
            end
            fprintf(']\n');
        elseif mod(k, 10) == 0
            fprintf('  Iteration %d: f = %.6f, |S| = %.6e, beta_avg = %.4f\n', ...
                k, f_k1, abs(S_k), mean(beta));
        end
        
        k = k + 1;
    end
    
    % Trim history
    alpha_history = alpha_history(:, 1:k-1);
    beta_history = beta_history(:, 1:k-1);
    f_history = f_history(1:k-1);
    g_history = g_history(1:k-1);
    S_history = S_history(1:k-1);
    
    % Return final values
    alpha_star = alpha;
    beta_star = beta;
    
    % Package history
    history.alpha = alpha_history;
    history.beta = beta_history;
    history.f = f_history;
    history.g = g_history;
    history.S = S_history;
end

% -------------------------------------------------------------------------
function alpha_new = update_alpha(alpha, beta, f, eta, x_o, a2)
% Update alpha using Nash best response equation (eq:alpha_k)
    n = length(alpha);
    alpha_new = zeros(n, 1);
    
    for i = 1:n
        % Best response: alpha_i = [(1+4a2*beta_i)/(2(1+2a2*beta_i))] * [1 - f + eta_i]
        coeff = (1 + 4*a2*beta(i)) / (2*(1 + 2*a2*beta(i)));
        alpha_new(i) = coeff * (1 - f + eta(i));
        
        % Project to [0,1]
        alpha_new(i) = max(0, min(1, alpha_new(i)));
    end
end

% -------------------------------------------------------------------------
function sum_val = consensus_sum(z_init, A_comm, epsilon, max_consensus_iter)
% Simulate consensus to compute sum (average * n)
% In distributed implementation, this would be actual message passing
    
    if nargin < 4
        max_consensus_iter = 100;
    end
    
    n = length(z_init);
    z = z_init;
    
    % Degree matrix
    D = diag(sum(A_comm, 2));
    
    for iter = 1:max_consensus_iter
        z_new = zeros(n, 1);
        for i = 1:n
            % Consensus update: z_i = z_i + epsilon * sum_j(z_j - z_i)
            neighbors = find(A_comm(i, :));
            z_new(i) = z(i) + epsilon * sum(z(neighbors) - z(i));
        end
        z = z_new;
        
        % Check convergence (all values equal)
        if max(abs(z - mean(z))) < 1e-8
            break;
        end
    end
    
    % Return sum (average * n)
    sum_val = mean(z) * n;
end

% -------------------------------------------------------------------------
function [alpha_nash, converged] = nash_equilibrium(beta, x_o, a2)
% Solve for Nash equilibrium given beta using equation (11) from paper
    
    n = length(x_o);
    x_r_o = sum(x_o);
    
    % Build matrix A and vector B from equation (11)
    A_mat = zeros(n, n);
    B_vec = zeros(n, 1);
    
    for i = 1:n
        for j = 1:n
            gamma_ij = x_o(j) / x_o(i);
            if i == j
                A_mat(i, j) = 2 * (1 + 2*a2*beta(i));
            else
                A_mat(i, j) = (1 + 4*a2*beta(i)) * gamma_ij;
            end
        end
        B_vec(i) = (1 + 4*a2*beta(i)) * sum(x_o) / x_o(i);
    end
    
    % Solve A * alpha = B
    alpha_nash = A_mat \ B_vec;
    
    % Check if converged (all alpha in [0,1])
    converged = all(alpha_nash >= 0 & alpha_nash <= 1);
    
    % Project to feasible set if needed
    alpha_nash = max(0, min(1, alpha_nash));
end

% -------------------------------------------------------------------------
function g = compute_g(alpha, eta, x_r_o)
% Compute aggregate profit metric g(alpha) = f*(1-f)*(x_r_o)^2
    f = sum(alpha .* eta);
    g = f * (1 - f) * x_r_o^2;
end

% -------------------------------------------------------------------------
function create_plots(algorithm_type, alpha_history, beta_history, ...
                     f_history, g_history, alpha_RN, g_RN, g_coop, n, eta, x_r_o)
% Create visualization plots

    if strcmp(algorithm_type, 'symmetric')
        % For symmetric case, just show final comparison
        figure('Position', [100, 100, 1000, 400]);
        
        % Comparison bar chart
        f_RN = sum(alpha_RN .* eta);
        strategies = {'Risk-Neutral', 'Pareto-Optimal', 'Cooperative'};
        g_values = [g_RN, g_history, g_coop];
        f_values = [f_RN, f_history, 0.5];
        
        subplot(1,2,1);
        bar(g_values);
        set(gca, 'XTickLabel', strategies);
        ylabel('Aggregate Profit g');
        title('Aggregate Profit Comparison');
        grid on;
        
        subplot(1,2,2);
        bar(f_values);
        hold on;
        yline(0.5, 'r--', 'LineWidth', 2);
        set(gca, 'XTickLabel', strategies);
        ylabel('Aggregate Fraction f');
        title('Bidding Fraction Comparison');
        legend('Actual', 'Optimal (0.5)');
        grid on;
        
    else
        % For iterative case, show convergence
        num_iter = length(f_history);
        
        figure('Position', [100, 100, 1400, 800]);
        
        % Plot 1: f convergence to 0.5
        subplot(2,3,1);
        plot(1:num_iter, f_history, 'b-', 'LineWidth', 2);
        hold on;
        yline(0.5, 'r--', 'LineWidth', 1.5);
        xlabel('Iteration');
        ylabel('f(\alpha)');
        title('Convergence of f to Pareto Condition');
        legend('f(\alpha)', 'Target f=0.5');
        grid on;
        
        % Plot 2: g convergence to maximum
        subplot(2,3,2);
        plot(1:num_iter, g_history, 'b-', 'LineWidth', 2);
        hold on;
        yline(g_coop, 'r--', 'LineWidth', 1.5);
        yline(g_RN, 'k:', 'LineWidth', 1.5);
        xlabel('Iteration');
        ylabel('g(\alpha)');
        title('Aggregate Profit Evolution');
        legend('g(\alpha)', 'Cooperative Max', 'Risk-Neutral Nash');
        grid on;
        
        % Plot 3: Alpha trajectories
        subplot(2,3,3);
        for i = 1:n
            plot(1:num_iter, alpha_history(i,:), 'LineWidth', 1.5);
            hold on;
        end
        yline(0.5, 'r--', 'LineWidth', 1);
        xlabel('Iteration');
        ylabel('\alpha_i');
        title('Bidding Fraction Trajectories');
        legend([arrayfun(@(i) sprintf('\\alpha_%d', i), 1:n, 'UniformOutput', false), 'Optimal']);
        grid on;
        
        % Plot 4: Beta trajectories
        subplot(2,3,4);
        for i = 1:n
            plot(1:num_iter, beta_history(i,:), 'LineWidth', 1.5);
            hold on;
        end
        xlabel('Iteration');
        ylabel('\beta_i');
        title('Risk Parameter Trajectories');
        legend(arrayfun(@(i) sprintf('\\beta_%d', i), 1:n, 'UniformOutput', false));
        grid on;
        
        % Plot 5: Error signal S(k)
        if isfield(struct('S', 0), 'S')
            subplot(2,3,5);
            semilogy(1:num_iter, abs(f_history - 0.5), 'b-', 'LineWidth', 2);
            xlabel('Iteration');
            ylabel('|S(k)| = |f - 0.5|');
            title('Error Signal (Log Scale)');
            grid on;
        end
        
        % Plot 6: Final comparison
        subplot(2,3,6);
        f_RN = sum(alpha_RN .* eta);
        strategies = {'Risk-Neutral', 'Pareto-Optimal', 'Cooperative'};
        g_values = [g_RN, g_history(end), g_coop];
        bar(g_values);
        set(gca, 'XTickLabel', strategies);
        ylabel('Aggregate Profit g');
        title('Final Comparison');
        text(2, g_history(end)*1.05, sprintf('%.1f%% improvement', ...
            100*(g_history(end)-g_RN)/g_RN), 'HorizontalAlignment', 'center');
        grid on;
    end
    
    sgtitle('Distributed Risk-Aware Nash Optimization Results', 'FontSize', 14, 'FontWeight', 'bold');
end

% -------------------------------------------------------------------------
function result = ternary(condition, true_val, false_val)
% Simple ternary operator for conditional string selection
    if condition
        result = true_val;
    else
        result = false_val;
    end
end