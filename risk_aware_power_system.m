
clc; clear;

%% parameters
mode = 2;   % 1 = centralized
            % 2 = decentralized
            % 3 = multiple runs for decentralized
n  = 3;
a2 = 0.8;
a1 = 1;
a0 = 0;

% Forecasted mean renewable generations x_i^o
x0 = [120; 100; 80];

kappa_l = 1.3; % load-to-renewable ratio
p_l = sum(x0)*kappa_l;

% Lower bound on beta: beta_i > -(4 a2)^(-1)
epsBeta = 1e-6;
beta_lb = -1/(4*a2) + epsBeta;
beta_prime = [0;0;-0.25];

rng(1);

%% ===================== RUN SELECTED MODE =====================
switch mode
    case 1
        run_centralized(n, a2, a1, a0, x0, beta_lb);

    case 2
        run_distributed(n, a2, a1, a0, x0, beta_lb, beta_prime);

    otherwise
        error('Unknown mode');
end

%% ===================== CENTRALIZED ROUTINE =====================
function run_centralized(n, a2, a1, a0, x0, beta_lb) %#ok<INUSD>
    fprintf('\n===== CENTRALIZED MODE =====\n');

    beta_RN = zeros(n,1);  % risk-neutral
    [alpha_RN, info_RN] = solve_alpha_nash(beta_RN, a2, x0);
    [g_RN, f_RN, xr0] = compute_g(alpha_RN, x0);

    fprintf('--- Risk-neutral baseline (beta=0) ---\n');
    disp(table((1:n)', beta_RN, alpha_RN, 'VariableNames', {'i','beta','alpha_star'}));
    fprintf('rcond(A)=%.3e, f=%.6f, g=%.6f (xr0=%.3f)\n\n', info_RN.rcondA, f_RN, g_RN, xr0);

    % Optimize beta to maximize g(beta) using fmincon on -g(beta)
    % Note: requires Optimization Toolbox.
    obj = @(beta) objective_neg_g(beta, a2, x0, beta_lb);

    beta0 = 0.01 * ones(n,1);
    lb    = beta_lb * ones(n,1);
    ub    = 10 * ones(n,1);

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

    fprintf('\n--- Optimized risk-aware betas (centralized) ---\n');
    disp(table((1:n)', beta_opt, alpha_opt, 'VariableNames', {'i','beta_opt','alpha_star_opt'}));
    fprintf('rcond(A)=%.3e, f=%.6f, g=%.6f\n', info_opt.rcondA, f_opt, g_opt2);

    % Optional 1D sweep plot for symmetric beta (beta_i all equal)
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
    plot(betas, gvals, 'LineWidth', 1.5); hold on;
    xlabel('\beta_i'); ylabel('g(\alpha^*(\beta))');
    title('g vs \beta_i');
    yline(g_RN, '--');
    scatter(beta_opt(1), g_opt, 'ro');
    legend("g", "RN baseline", "g(\beta^*)", 'Location','best');
    grid on;

    figure;
    plot(betas, fvals, 'LineWidth', 1.5); hold on;
    xlabel('\beta_i'); ylabel('f(\alpha^*(\beta))');
    title('f vs \beta_i');
    yline(0.5, '--');
    scatter(beta_opt(1), f_opt, 'ro');
    legend("f", "f(\alpha^*)=0.5", "f(\alpha^*(\beta^*))", 'Location','best');
    grid on;
end

%% ===================== DISTRIBUTED ROUTINE =====================
function [beta_path, alpha_path] = run_distributed(n, a2, a1, a0, x0, beta_lb, beta_prime)
    fprintf('\n===== DISTRIBUTED MODE =====\n');

    % Inner loop (solve alpha* for fixed beta_prime)
    K_alpha = 2000;
    tau     = 0.2;
    tol_a   = 1e-10;

    % Outer loop (adjust beta_prime to enforce Theorem 2 condition f(alpha*)=1/2)
    K_beta = 200;
    kappa  = 0.6;
    tol_f  = 1e-6;

    % Constraints on beta_prime: beta_prime >= -1/4 corresponds to beta >= -(4a2)^(-1)
    beta_prime_min = -4*a2^-1;
    beta_prime_max =  10;

    xr0  = sum(x0);
    eta  = x0 / xr0;

    proj01 = @(z) min(1, max(0, z));
    projB  = @(z) min(beta_prime_max, max(beta_prime_min, z));

    hist.f = zeros(K_beta,1);
    hist.g = zeros(K_beta,1);
    hist.sum_beta_prime = zeros(K_beta,1);
    beta_hist = zeros(K_beta, n);
    alpha = 0.5*ones(n,1);
    alpha_star = 0.5*ones(n,1);

    for t = 1:K_beta
        f = eta' * alpha;
        alpha_map = ((1 + 4*beta_prime) ./ eta) * (1 - f);
        alpha_next = proj01((1 - tau)*alpha + tau*alpha_map);
        alpha = alpha_next;

        alpha_star = alpha;
        f_star = eta' * alpha_star;
        [g_star, ~, ~] = compute_g(alpha_star, x0);

        hist.f(t) = f_star;
        hist.g(t) = g_star;
        hist.sum_beta_prime(t) = sum(beta_prime);
        beta_hist(t,:) = beta_prime.';

        % stop if Theorem 2 condition achieved
        if abs(f_star - 0.5) < tol_f
            fprintf('Converged at outer iter %d: f = %.6f\n', t, f_star);
            hist.f = hist.f(1:t);
            hist.g = hist.g(1:t);
            hist.sum_beta_prime = hist.sum_beta_prime(1:t);
            beta_hist = beta_hist(1:t,:);
            break;
        end

        % Estimate S from f: f = S/(1+S) => S = f/(1-f)
        if f_star >= 1 - 1e-12
            S_hat = 1e12;
        else
            S_hat = f_star / (1 - f_star);
        end

        % outer update: drive S -> 1  (equivalently sum beta_prime -> (1-n)/4)
        delta = (1 - S_hat) / (4*n);
        beta_prime = projB(beta_prime + kappa * delta * ones(n,1));
    end

    % Centralized check for the same beta_prime:
    % Convert to beta for A\B build: beta_prime = a2 * beta => beta = beta_prime / a2
    beta = beta_prime / a2;
    [alpha_central, info] = solve_alpha_nash(beta, a2, x0);
    f_central = (x0/sum(x0))' * alpha_central;

    fprintf('\nFinal results (distributed):\n');
    fprintf('beta_prime = [%s]^T\n', num2str(beta_prime','%.6f '));
    fprintf('sum beta_prime = %.6f (target = (1-n)/4 = %.6f)\n', sum(beta_prime), (1-n)/4);

    fprintf('alpha_star (distributed) = [%s]^T\n', num2str(alpha_star','%.6f '));
    fprintf('f(alpha_star) = %.6f (target = 0.5)\n', eta'*alpha_star);
    fprintf('g(alpha_star) = %.6f', compute_g(alpha_star, x0))

    % Plots
    figure; plot(hist.f,'LineWidth',1.5); grid on;
    xlabel('Iterations'); ylabel('f(\alpha^*)');
    title('f(\alpha^*) over Iterationss');
    yline(0.5,'--');

    figure; plot(hist.g,'LineWidth',1.5); grid on;
    xlabel('Iterations'); ylabel('g(\alpha^*)');
    title('g(\alpha^*) over Iterationss');

    figure; plot(hist.sum_beta_prime,'LineWidth',1.5); grid on;
    xlabel('Iterations'); ylabel('\Sigma_i \beta_i''');
    title('\Sigma_i \beta_i'''' over Iterationss');
    yline((1-n)/4,'--');

    figure; plot(beta_hist,'LineWidth',1.5); grid on;
    xlabel('Iterations'); ylabel('\beta_i''');
    title('\beta_i'''' trajectories');
    legend(arrayfun(@(i) sprintf('\\beta_{%d}''''', i), 1:n, 'UniformOutput', false), ...
           'Location','best');
end

%% ===================== COMMON HELPER FUNCTIONS =====================
function [alpha_star, info] = solve_alpha_nash(beta, a2, x0)
    % Solves A(beta)*alpha = B(beta) from Theorem 1 / Eq (11) (general n)
    n = numel(beta);
    gamma = x0' ./ x0; % gamma(i,j) = x0(j)/x0(i)

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
    xr0 = sum(x0);
    eta = x0 / xr0;
    f = sum(alpha(:) .* eta(:));
    g = f*(1 - f) * (xr0^2);
end

function val = objective_neg_g(beta, a2, x0, beta_lb)
    [alpha, ok] = safe_solve_alpha(beta, a2, x0, beta_lb);
    if ~ok
        val = 1e30;
        return;
    end
    [g, ~, ~] = compute_g(alpha, x0);
    val = -g;
end

function [alpha, ok] = safe_solve_alpha(beta, a2, x0, beta_lb)
    ok = true;

    if any(beta <= beta_lb)
        ok = false;
        alpha = nan(size(beta));
        return;
    end

    [alpha, info] = solve_alpha_nash(beta, a2, x0);

    if any(~isfinite(alpha)) || info.rcondA < 1e-12
        ok = false;
        alpha = nan(size(beta));
        return;
    end
end
