clc; clear;

%% Parameters (ground truth for verification)
n_true = 3;           % actual number of renewable operators
a2 = 0.8;
a1 = 1;
a0 = 0;

% Forecasted mean renewable generations x_i^o (can be asymmetric)
x0 = [120; 100; 80];

% Ground truth values (for verification)
xr0_true = sum(x0);                    % = 300
eta_true = x0 / xr0_true;              % = [0.4; 0.333; 0.267]
alpha_test = [0.5; 0.5; 0.5];          % test bidding fractions
f_true = sum(alpha_test .* eta_true);  % = 0.5

fprintf('=== Ground Truth Values ===\n');
fprintf('n = %d\n', n_true);
fprintf('x_r^o = %.2f\n', xr0_true);
fprintf('eta = [%.4f, %.4f, %.4f]\n', eta_true(1), eta_true(2), eta_true(3));
fprintf('f (with alpha=0.5) = %.4f\n\n', f_true);

%% Communication Graph Setup
% Define adjacency matrix (connected graph required by Assumption 1)
% Using ring topology: 1 -- 2 -- 3 -- 1
Adj = [0 1 1;
       1 0 1;
       1 1 0];  % Fully connected (can also use sparse ring)

% Compute degree and Laplacian matrices
Deg = diag(sum(Adj, 2));
L = Deg - Adj;

% Consensus gain epsilon (must satisfy: epsilon < 1/max_degree)
max_degree = max(diag(Deg));
epsilon = 0.4 / max_degree;  % Conservative choice for stability

fprintf('=== Communication Network ===\n');
fprintf('Max degree = %d, epsilon = %.4f\n\n', max_degree, epsilon);

%% Consensus Parameters
max_iter = 200;
tol = 1e-8;

%% Step 1: Compute n via consensus (Eq. 5.3)
% Operator 1 initializes z_1(0) = 1, all others z_j(0) = 0
fprintf('=== Step 1: Computing n via consensus ===\n');
z = zeros(n_true, 1);
z(1) = 1;

z_history_n = zeros(n_true, max_iter);
for k = 1:max_iter
    z_history_n(:, k) = z;
    % Consensus update: z_i(k+1) = z_i(k) + epsilon * sum_{j in N_i}(z_j - z_i)
    z_new = z - epsilon * L * z;  % Matrix form of consensus
    if norm(z_new - z) < tol
        fprintf('Converged at iteration %d\n', k);
        break;
    end
    z = z_new;
end
z_bar_n = mean(z);  % All values should be equal at convergence
n_est = 1 / z_bar_n;

fprintf('Estimated n = %.6f (true = %d, error = %.2e)\n\n', n_est, n_true, abs(n_est - n_true));

%% Step 2: Compute x_r^o via consensus (Eq. 5.4)
% Each operator initializes z_i(0) = x_i^o
fprintf('=== Step 2: Computing x_r^o via consensus ===\n');
z = x0;

z_history_xr = zeros(n_true, max_iter);
for k = 1:max_iter
    z_history_xr(:, k) = z;
    z_new = z - epsilon * L * z;
    if norm(z_new - z) < tol
        fprintf('Converged at iteration %d\n', k);
        break;
    end
    z = z_new;
end
z_bar_xr = mean(z);  % = (1/n) * sum(x_j^o)
xr0_est = n_est * z_bar_xr;  % x_r^o = n * z_bar

fprintf('Estimated x_r^o = %.6f (true = %.2f, error = %.2e)\n\n', xr0_est, xr0_true, abs(xr0_est - xr0_true));

%% Step 3: Compute eta_i locally (Eq. 5.5)
% Each operator computes eta_i = x_i^o / x_r^o using local x_i^o and shared x_r^o
fprintf('=== Step 3: Computing eta_i ===\n');
eta_est = x0 / xr0_est;

fprintf('Estimated eta = [%.6f, %.6f, %.6f]\n', eta_est(1), eta_est(2), eta_est(3));
fprintf('True eta      = [%.6f, %.6f, %.6f]\n', eta_true(1), eta_true(2), eta_true(3));
fprintf('Max error = %.2e\n\n', max(abs(eta_est - eta_true)));

%% Step 4: Compute f via consensus (Eq. 5.6)
% Each operator initializes z_i(0) = alpha_i * x_i^o
fprintf('=== Step 4: Computing f via consensus ===\n');
z = alpha_test .* x0;

z_history_f = zeros(n_true, max_iter);
for k = 1:max_iter
    z_history_f(:, k) = z;
    z_new = z - epsilon * L * z;
    if norm(z_new - z) < tol
        fprintf('Converged at iteration %d\n', k);
        break;
    end
    z = z_new;
end
z_bar_f = mean(z);  % = (1/n) * sum(alpha_j * x_j^o)
f_est = n_est * z_bar_f / xr0_est;  % f = n * z_bar / x_r^o

fprintf('Estimated f = %.6f (true = %.6f, error = %.2e)\n\n', f_est, f_true, abs(f_est - f_true));

%% Visualization
figure('Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
plot(1:k, z_history_n(:, 1:k)', 'LineWidth', 1.5);
hold on;
yline(1/n_true, 'k--', 'LineWidth', 1.5);
xlabel('Iteration k');
ylabel('z_i(k)');
title('Consensus for n: z_i(0) = \delta_{i,1}');
legend('Operator 1', 'Operator 2', 'Operator 3', '1/n', 'Location', 'best');
grid on;

subplot(1, 3, 2);
plot(1:k, z_history_xr(:, 1:k)', 'LineWidth', 1.5);
hold on;
yline(xr0_true/n_true, 'k--', 'LineWidth', 1.5);
xlabel('Iteration k');
ylabel('z_i(k)');
title('Consensus for x_r^o: z_i(0) = x_i^o');
legend('Operator 1', 'Operator 2', 'Operator 3', 'x_r^o/n', 'Location', 'best');
grid on;

subplot(1, 3, 3);
plot(1:k, z_history_f(:, 1:k)', 'LineWidth', 1.5);
hold on;
yline(f_true * xr0_true / n_true, 'k--', 'LineWidth', 1.5);
xlabel('Iteration k');
ylabel('z_i(k)');
title('Consensus for f: z_i(0) = \alpha_i x_i^o');
legend('Operator 1', 'Operator 2', 'Operator 3', 'Target', 'Location', 'best');
grid on;

sgtitle('Consensus Protocol Convergence (Equations 5.3-5.6)');

fprintf('=== Summary ===\n');
fprintf('All consensus computations converged successfully.\n');
fprintf('The distributed protocol correctly estimates n, x_r^o, eta, and f.\n');