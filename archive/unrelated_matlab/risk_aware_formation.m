%% risk_beta_formation.m
% Beta-Parameterized Risk-Aware Distributed Formation Control
%
% USAGE:
%   Set beta at the top, then run.
%
%   beta = -1  risk-seeking  (high gains, fast transient, more variance)
%   beta =  0  risk-neutral  (standard consensus, k = k_base always)
%   beta = +1  risk-averse   (low gains, slow transient, less variance)
%
% How beta works:
%   k_i(t) = k_base / (1 + beta * sigma_i(t)^2)
%
%   sigma_i is a running EMA of each agent's CONSENSUS error —
%   how far agent i is from where its neighbors expect it to be.
%   This goes to zero at steady state for all beta, so the gain
%   converges to k_base. Beta only shapes the transient behavior.
%
% UAV 3 scenario:
%   t < t_fail    : all 4 UAVs hold diamond
%   t_fail→rejoin : UAV 3 loses comms, drifts freely (pure noise, no u)
%                   UAVs 1,2,4 hold their original diamond slots
%   t >= t_rejoin : UAV 3 regains comms from wherever it drifted,
%                   consensus pulls it back into the diamond

clear; clc; close all;

%% ── SET BETA HERE ────────────────────────────────────────────────────────
beta = -1;   % try: -1.0, -0.5, 0.0, 0.5, 1.0
%──────────────────────────────────────────────────────────────────────────

%% ── Parameters ───────────────────────────────────────────────────────────
dt            = 0.05;
T_end         = 30;
t             = 0:dt:T_end;
N_steps       = length(t);
n_UAV         = 4;

k_base        = 1.5;
k_min         = 0.3;    % tighter floor — prevents runaway for beta<0
k_max         = 4.0;    % tighter ceiling than before

sigma_process = 0.10;   % process noise std dev [m] during normal flight
sigma_drift   = 0.40;   % UAV 3 drift noise while comms are down [m]

% EMA smoothing for sigma estimate
% sigma_i^2(k) = (1-alpha)*sigma_i^2(k-1) + alpha*e_consensus_i^2(k)
ema_alpha     = 0.05;   % smaller = longer memory, smoother sigma

t_fail        = 5.0;
t_rejoin      = 14.0;

%% ── Formation geometry ───────────────────────────────────────────────────
formation_pos = [ 0,  2;   % UAV 1  top
                 -2,  0;   % UAV 2  left
                  0, -2;   % UAV 3  bottom
                  2,  0];  % UAV 4  right

d = zeros(n_UAV, n_UAV, 2);
for i = 1:n_UAV
    for j = 1:n_UAV
        d(i,j,:) = formation_pos(i,:) - formation_pos(j,:);
    end
end

idx3 = [1, 2, 4];   % active UAVs during dead window

%% ── Storage ──────────────────────────────────────────────────────────────
X_hist      = zeros(N_steps, n_UAV, 2);
Eerr_hist   = zeros(N_steps, 1);
K_hist      = zeros(N_steps, n_UAV);
Sigma_hist  = zeros(N_steps, n_UAV);
Jtot_hist   = zeros(N_steps, 1);

%% ── Reproducible noise ───────────────────────────────────────────────────
rng(42);
base_noise = randn(N_steps, n_UAV, 2);

%% ── Initial conditions ───────────────────────────────────────────────────
rng(7);
directions = rand(n_UAV,2) - 0.5;
directions = directions ./ vecnorm(directions, 2, 2);  % unit vectors
distance   = 3.0;                                       % meters from slot
x = formation_pos + distance * directions;

% Sigma initialised from starting consensus error
Sigma = zeros(n_UAV, 1);
for i = 1:n_UAV
    e_consensus = 0;
    for j = 1:n_UAV
        if i ~= j
            des         = squeeze(d(i,j,:))';
            e_consensus = e_consensus + norm((x(i,:) - x(j,:)) - des);
        end
    end
    Sigma(i) = e_consensus / (n_UAV - 1);
end

uav3_dead = false;
uav3_back = false;

%% ── Simulation loop ──────────────────────────────────────────────────────
fprintf('Simulating  beta = %+.2f\n', beta);

for k = 1:N_steps
    tk = t(k);

    %% ── Events ───────────────────────────────────────────────────────────
    if tk >= t_fail && ~uav3_dead
        uav3_dead = true;
        fprintf('  [t=%5.1fs]  UAV 3 comms LOST\n', tk);
    end

    if tk >= t_rejoin && uav3_dead && ~uav3_back
        uav3_back = true;
        uav3_dead = false;
        drift_dist = norm(x(3,:) - formation_pos(3,:));
        % Sigma_3 at rejoin = mean consensus error to its neighbors
        e_rj = 0;
        for j = idx3
            des  = squeeze(d(3,j,:))';
            e_rj = e_rj + norm((x(3,:) - x(j,:)) - des);
        end
        Sigma(3) = e_rj / length(idx3);
        fprintf('  [t=%5.1fs]  UAV 3 comms BACK — %.2fm from slot, sigma=%.2f\n', ...
                tk, drift_dist, Sigma(3));
    end

    %% ── Gains ────────────────────────────────────────────────────────────
    K = compute_beta_gains(Sigma, beta, k_base, k_min, k_max, n_UAV);

    %% ── Control inputs ───────────────────────────────────────────────────
    u = zeros(n_UAV, 2);

    if uav3_dead
        % UAVs 1,2,4 hold original diamond slots
        for ii = 1:length(idx3)
            i   = idx3(ii);
            u_i = zeros(1,2);
            for jj = 1:length(idx3)
                if ii ~= jj
                    j    = idx3(jj);
                    k_ij = (K(i) + K(j)) / 2;
                    des  = squeeze(d(i,j,:))';
                    u_i  = u_i - k_ij * (x(i,:) - x(j,:) - des);
                end
            end
            u(i,:) = u_i;
        end
        % UAV 3: u = 0, drifts freely

    else
        for i = 1:n_UAV
            u_i = zeros(1,2);
            for j = 1:n_UAV
                if i ~= j
                    k_ij = (K(i) + K(j)) / 2;
                    des  = squeeze(d(i,j,:))';
                    u_i  = u_i - k_ij * (x(i,:) - x(j,:) - des);
                end
            end
            u(i,:) = u_i;
        end
    end

    %% ── Integrate ────────────────────────────────────────────────────────
    w = zeros(n_UAV, 2);
    for i = 1:n_UAV
        if uav3_dead && i == 3
            w(i,:) = sigma_drift   * squeeze(base_noise(k,i,:))';
        else
            w(i,:) = sigma_process * squeeze(base_noise(k,i,:))';
        end
    end
    x = x + dt*(u + w);

    %% ── Update sigma via CONSENSUS error (not absolute position error) ───
    % Consensus error for agent i = average mismatch to each neighbor
    % relative to desired offset. This → 0 when formation is achieved,
    % regardless of noise floor. Safe for beta < 0.
    for i = 1:n_UAV
        if uav3_dead && i == 3
            % UAV 3 is offline — freeze its sigma during dead window
            continue
        end

        active_neighbors = setdiff(1:n_UAV, i);
        if uav3_dead
            active_neighbors = setdiff(idx3, i);
        end

        e_sum = 0;
        for j = active_neighbors
            des   = squeeze(d(i,j,:))';
            e_sum = e_sum + norm((x(i,:) - x(j,:)) - des);
        end
        e_consensus = e_sum / length(active_neighbors);

        % EMA update
        Sigma(i) = sqrt((1-ema_alpha)*Sigma(i)^2 + ema_alpha*e_consensus^2);
    end

    %% ── Metrics ──────────────────────────────────────────────────────────
    active = 1:n_UAV;
    if uav3_dead, active = idx3; end

    errs = zeros(length(active),1);
    for ii = 1:length(active)
        i = active(ii);
        errs(ii) = norm(x(i,:) - formation_pos(i,:));
    end

    mean_err = mean(errs);
    var_err  = var(errs);

    %% ── Store ────────────────────────────────────────────────────────────
    X_hist(k,:,:)   = x;
    Eerr_hist(k)    = mean_err;
    K_hist(k,:)     = K';
    Sigma_hist(k,:) = Sigma';
    Jtot_hist(k)    = mean_err^2 + beta*var_err;
end

fprintf('Done.\n');

%% ── Plot ─────────────────────────────────────────────────────────────────
plot_beta_results(t, dt, beta, X_hist, Eerr_hist, K_hist, Sigma_hist, ...
    Jtot_hist, formation_pos, n_UAV, t_fail, t_rejoin);