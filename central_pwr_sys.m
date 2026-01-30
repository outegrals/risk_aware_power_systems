
clc; clear;

%% parameters
n  = 3;                % number of renewable operators
a2 = 0.8;
a1 = 1;
a0 = 0;

% Forecasted mean renewable generations x_i^o (can be asymmetric)
x0 = [120; 100; 80];

epsBeta = 1e-6;
beta_lb = -1/(4*a2) + epsBeta;
rng(1);


beta_RN = zeros(n,1);  % risk-neutral
[alpha_RN, info_RN] = solve_alpha_nash(beta_RN, a2, x0);
[g_RN, f_RN, xr0] = compute_g(alpha_RN, x0);
fprintf('--- Risk-neutral baseline (beta=0) ---\n');
disp(table((1:n)', beta_RN, alpha_RN, 'VariableNames', {'i','beta','alpha_star'}));
fprintf('f = %.6f, g = %.6f (xr0=%.3f)\n\n', f_RN, g_RN, xr0);

%% beta optimization to maximize g(beta)

%todo: use analytical approach if possible?
% easy to do with n=2, but n>2 is more complex

% fmincon: minimize negative g(beta)
obj = @(beta) objective_neg_g(beta, a2, x0, beta_lb);
beta0 = 0.01 * ones(n,1);           % initial guess
lb    = beta_lb * ones(n,1);
ub    = 10 * ones(n,1);             % arbitrary large cap; adjust as desired
opts = optimoptions('fmincon', ...
    'Display','iter', ...
    'Algorithm','sqp', ...
    'MaxFunctionEvaluations', 5e4, ...
    'OptimalityTolerance', 1e-10, ...
    'StepTolerance', 1e-12);
[beta_opt, fval] = fmincon(obj, beta0, [], [], [], [], lb, ub, [], opts);
g_opt = -fval;

[alpha_opt, info_opt] = solve_alpha_nash(beta_opt, a2, x0);
[g_opt2, f_opt, ~] = compute_g(alpha_opt, x0);
fprintf('\n--- Optimized risk-aware betas ---\n');
disp(table((1:n)', beta_opt, alpha_opt, 'VariableNames', {'i','beta_opt','alpha_star_opt'}));
fprintf('f = %.6f, g = %.6f\n', f_opt, g_opt2);

%% plotting
betas = linspace(beta_lb, 2, 400);
gvals = nan(size(betas));
fvals = nan(size(betas));
for k = 1:numel(betas)
    beta = betas(k) * ones(n,1);
    [alpha, ok] = safe_solve_alpha(beta, a2, x0, beta_lb);
    if ~ok, continue; end
    [g, f, ~] = compute_g(alpha, x0);
    gvals(k) = g; fvals(k) = f;
end
figure; 
plot(betas, gvals, 'LineWidth', 1.5);
hold on;
xlabel('\beta_i'); ylabel('g(\alpha^*(\beta))');
title('g vs \beta_i');
yline(g_RN, '--');
scatter( beta_opt(1), g_opt, 'ro');
legend("g", "RN Baseline", "g(\beta*)");
grid on;

figure; 
plot(betas, fvals, 'LineWidth', 1.5);
hold on;
xlabel('\beta_i'); ylabel('f(\alpha^*(\beta))');
title('f vs \beta_i');
yline(0.5, '--'); grid on;
scatter( beta_opt(1), f_opt, 'ro');
legend("f", "f(\alpha*)=0.5", "f(\alpha*(\beta*))");

%% helper functions
function [alpha_star, info] = solve_alpha_nash(beta, a2, x0)
    % Solves A(beta)*alpha = B(beta) from Theorem 1 / Eq (11) (general n)
    % beta: nx1, x0: nx1, a2 scalar > 0

    n = numel(beta);
    gamma = zeros(n,n);

    for i = 1:n
        for j = 1:n
            gamma(i,j) = x0(j)/x0(i);  % gamma_ij = x_j^o / x_i^o
        end
    end

    A = zeros(n,n);
    B = zeros(n,1);

    for i = 1:n
        A(i,i) = 2*(1 + 2*a2*beta(i));
        for j = 1:n
            if j ~= i
                A(i,j) = (1 + 4*a2*beta(i)) * gamma(i,j);
            end
        end
        B(i) = (1 + 4*a2*beta(i)) * sum(gamma(i,:));
    end

    alpha_star = A \ B;

    info.A = A;
    info.B = B;
    info.rcondA = rcond(A);
end

function [g, f, xr0] = compute_g(alpha, x0)
    % Implements Eqs (12)-(13):
    % xr0 = sum x_i^o
    % eta_i = x_i^o / xr0
    % f = sum alpha_i * eta_i
    % g = f(1-f) (xr0)^2

    xr0 = sum(x0);
    eta = x0 / xr0;

    f = sum(alpha(:) .* eta(:));
    g = f*(1 - f) * (xr0^2);
end

function val = objective_neg_g(beta, a2, x0, beta_lb)
    % fmincon objective: minimize -g(beta)
    [alpha, ok] = safe_solve_alpha(beta, a2, x0, beta_lb);
    if ~ok
        val = 1e30;  % penalize infeasible / ill-conditioned
        return;
    end
    [g, ~, ~] = compute_g(alpha, x0);
    val = -g;
end

function [alpha, ok] = safe_solve_alpha(beta, a2, x0, beta_lb)
    % Enforce the admissible range and numerical stability.
    ok = true;

    if any(beta <= beta_lb)
        ok = false;
        alpha = nan(size(beta));
        return;
    end

    [alpha, info] = solve_alpha_nash(beta, a2, x0);

    % Basic numerical sanity checks
    if any(~isfinite(alpha)) || info.rcondA < 1e-12
        ok = false;
        alpha = nan(size(beta));
        return;
    end
end
