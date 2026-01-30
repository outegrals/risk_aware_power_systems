%% test_dist_theorem2.m
% Distributed risk-aware Nash: Lemma 2 (inner alpha solve) + Theorem 2 (outer beta optimize)
% Also compares against centralized alpha* = A\B for the same beta'.

clear; clc;

%% ---------------- User settings ----------------
n  = 3;                         % number of renewable operators
a2 = 0.8;
a1 = 1;
a0 = 0;                     % cost coefficient (only used to map beta' <-> beta)
x0 = [120; 100; 80];              % expected renewable outputs x_i^o  (n x 1)
% NOTE: only x0 matters for alpha/beta updates here, not Sigma.

% Initialization for beta (risk attitude)
%beta_prime = zeros(n,1);        % beta'_i = a2 * beta_i
beta_prime = [-1; -1; -1];
% Example: start slightly risk-neutral, or slightly risk-seeking:
% beta_prime = -0.01*ones(n,1);

% Inner loop settings (solve alpha* for fixed beta')
K_alpha = 2000;     % max inner iterations
tau     = 0.2;      % relaxation for alpha update (0<tau<=1)
tol_a   = 1e-10;    % inner loop stopping tolerance

% Outer loop settings (optimize beta' to maximize profit => f(alpha*) = 1/2)
K_beta = 200;       % max Iterations
kappa  = 0.6;       % stepsize for beta' update (0<kappa<2 recommended)
tol_f  = 1e-6;      % stopping tolerance for |f - 1/2|

% Constraints (optional but recommended)
beta_prime_min = -0.24;   % corresponds to beta_i >= -(4 a2)^{-1} (i.e., beta'_i >= -1/4)
beta_prime_max =  10;     % arbitrary cap

%% ---------------- Precompute constants ----------------
xr0  = sum(x0);
eta  = x0 / xr0;                % eta_i = x_i^o / x_r^o

% Helper: projection to [0,1]
proj01 = @(z) min(1, max(0, z));
% Helper: projection for beta'
projB = @(z) min(beta_prime_max, max(beta_prime_min, z));

%% ---------------- Outer loop: update beta' ----------------
alpha_star = zeros(n,1);

hist.f = zeros(K_beta,1);
hist.S = zeros(K_beta,1);
hist.g = zeros(K_beta,1);
hist.sum_beta_prime = zeros(K_beta,1);
alpha = 0.5*ones(n,1); % initial alpha for inner loop (any in [0,1] works)
beta_hist = zeros(K_beta, n);   % store beta'_i over time
for t = 1:K_beta
    % ---- Inner loop: compute alpha*(beta') via Lemma 2 fixed-point iteration ----
    
    %for k = 1:K_alpha
        f = eta' * alpha;  % scalar aggregate f(alpha)
        alpha_map = ((1 + 4*beta_prime) ./ eta) * (1 - f); % componentwise map (may go out of [0,1])

        alpha_next = proj01((1 - tau)*alpha + tau*alpha_map);

        alpha = alpha_next;
    %end

    alpha_star = alpha;
    f_star = eta' * alpha_star;
    [g, ~, ~] = compute_g(alpha_star, x0);

    % Record
    hist.f(t) = f_star;
    hist.g(t) = g;

    % Convert f -> S estimate using f = S/(1+S) => S = f/(1-f)
    % Guard against division by zero if f ~ 1
    if f_star >= 1 - 1e-12
        S_hat = 1e12;
    else
        S_hat = f_star / (1 - f_star);
    end
    hist.S(t) = S_hat;
    hist.sum_beta_prime(t) = sum(beta_prime);

    % Stop if Theorem 2 condition is met (f = 1/2)
    if abs(f_star - 0.5) < tol_f
        fprintf('Converged at outer iter %d: f = %.6f\n', t, f_star);
        hist.f = hist.f(1:t);
        hist.S = hist.S(1:t);
        hist.g = hist.g(1:t);
        hist.sum_beta_prime = hist.sum_beta_prime(1:t);
        beta_hist = beta_hist(1:t-1, :);
        break;
    end

    % ---- Outer update: drive S -> 1 (equivalently f -> 1/2) ----
    % Target S* = 1  <=>  sum beta'_i = (1-n)/4
    % Simple fair distributed rule: everyone shifts equally based on scalar error (1 - S_hat)
    delta = (1 - S_hat) / (4*n);           % equal share so S_{t+1} = S_t + kappa(1-S_t)
    beta_prime = projB(beta_prime + kappa * delta * ones(n,1));
    %beta_bar   = mean(beta_prime);
    %beta_prime = beta_bar*ones(n,1);
    beta_hist(t,:) = beta_prime.';
end
%beta_hist = beta_hist(1:t-1,:);

%% ---------------- Centralized check for final beta' ----------------
% Build A and B (Theorem 1 form with beta' = a2 beta):
% gamma_ij = x_j^o / x_i^o
gamma = x0' ./ x0; % produces n x n matrix with gamma(i,j)=x0(j)/x0(i)

A = zeros(n,n);
B = zeros(n,1);
for i = 1:n
    for j = 1:n
        if i == j
            A(i,j) = 2*(1 + 2*beta_prime(i));
        else
            A(i,j) = (1 + 4*beta_prime(i)) * gamma(i,j);
        end
    end
    B(i) = (1 + 4*beta_prime(i)) * sum(gamma(i,:));
end

alpha_central = A \ B;

%% ---------------- Report ----------------
fprintf('\nFinal results:\n');
fprintf('beta_prime (distributed) = [%s]^T\n', num2str(beta_prime','%.6f '));
fprintf('sum beta_prime = %.6f (target = (1-n)/4 = %.6f)\n', sum(beta_prime), (1-n)/4);

fprintf('alpha_star = [%s]^T\n', num2str(alpha_star','%.6f '));
fprintf('f(alpha_star) = %.6f \n', eta'*alpha_star);
fprintf('g(alpha_star) = %.6f', compute_g(alpha_star, x0))

%% ---------------- Optional plots ----------------
figure; plot(hist.f,'LineWidth',1.5); grid on;
xlabel('Iteration'); ylabel('f(\alpha^*)');
title('f(\alpha^*) iterations');


figure; plot(hist.g,'LineWidth',1.5); grid on;
xlabel('Iteration'); ylabel('g');
title('g(\alpha^*) iterations');

figure; plot(hist.sum_beta_prime,'LineWidth',1.5); grid on;
xlabel('Iteration'); ylabel('\Sigma_i \beta_i''');
title('Sum of beta iterations');

figure;
plot(beta_hist, 'LineWidth', 1.5);
grid on;
xlabel('Iteration');
ylabel('\beta_i''');
title('\beta_i iterations');

legend(arrayfun(@(i) sprintf('\\beta_{%d}''', i), 1:n, ...
       'UniformOutput', false), 'Location', 'best');

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