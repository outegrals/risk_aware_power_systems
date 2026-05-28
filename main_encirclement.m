%% main_encirclement.m
% Risk-Aware Distributed UAV Encirclement with Cooperative Recovery
%
% Three cases:
%   Case 1 — No comms loss,  fixed beta
%   Case 2 — UAV3 comms loss, fixed beta
%   Case 3 — UAV3 comms loss, online beta (gradient descent)

clear; clc; close all;

%% ══════════════════════════════════════════════════════════════════════════
%  SET TARGET MODE HERE
%  'lissajous' — smooth figure-8 (deterministic, good for paper figures)
%  'random'    — smooth random walk (more realistic, stresses the controller)
% ══════════════════════════════════════════════════════════════════════════
p.target_mode = 'lissajous';   % <── change me

%% ── Lissajous (figure-8) params ──────────────────────────────────────────
p.target_Ax = 4.0;    % x amplitude [m]
p.target_Ay = 2.5;    % y amplitude [m]
p.target_wx = 0.20;   % x frequency [rad/s]
p.target_wy = 0.40;   % y frequency [rad/s]  (2:1 ratio → figure-8)

%% ── Random target params ─────────────────────────────────────────────────
p.target_speed  = 1.2;  % max speed [m/s]
p.target_tau    = 3.0;  % smoothing time constant [s]  (larger = lazier turns)
p.target_bounds = 7.0;  % soft boundary radius [m]
p.seed_target   = 15;   % fixed seed — all 3 cases share same path

%% ── Fixed beta (Cases 1 and 2) ───────────────────────────────────────────
beta_fixed = 0.0;

%% ── Online beta params (Case 3) ─────────────────────────────────────────
beta_init = 0.0;
eta       = 0.15;

%% ── Shared simulation params ─────────────────────────────────────────────
p.dt      = 0.05;
p.T_end   = 30;
p.t       = 0 : p.dt : p.T_end;
p.N_steps = length(p.t);
p.n_UAV   = 4;

p.r0           = 3.0;
p.k_radial     = 1.8;
p.k_theta      = 1.2;
p.k_track      = 1.0;

p.sigma_process = 0.08;
p.ema_alpha     = 0.06;

p.k_base  = 1.5;
p.k_min   = 0.3;
p.k_max   = 4.0;

p.t_fail   = 8.0;
p.t_rejoin = 18.0;

p.seed_noise = 42;
p.seed_ic    = 7;

%% ── Run ──────────────────────────────────────────────────────────────────
fprintf('Target mode: %s\n\n', p.target_mode);

fprintf('=== Case 1: No comms loss  (beta = %+.2f) ===\n', beta_fixed);
p.beta   = beta_fixed;
p.online = false;
results1 = run_encirclement(p, false);

fprintf('\n=== Case 2: No comms loss  (online beta, init=%.2f, eta=%.3f) ===\n', ...
        beta_init, eta);
p.beta   = beta_init;
p.eta    = eta;
p.online = true;
results2 = run_encirclement(p, false);

fprintf('\n=== Case 3: UAV3 comms loss  (beta = %+.2f) ===\n', beta_fixed);
p.beta   = beta_fixed;
p.online = false;
results3 = run_encirclement(p, true);

fprintf('\n=== Case 4: UAV3 comms loss  (online beta, init=%.2f, eta=%.3f) ===\n', ...
        beta_init, eta);
p.beta   = beta_init;
p.eta    = eta;
p.online = true;
results4 = run_encirclement(p, true);

%% ── Plot ─────────────────────────────────────────────────────────────────
%plot_encirclement(p, results1, results2, results3);