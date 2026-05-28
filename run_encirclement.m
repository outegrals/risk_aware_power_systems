function results = run_encirclement(p, comm_loss)
% RUN_ENCIRCLEMENT  Distributed UAV encirclement of a smooth random target.
%
% TARGET: pre-generated smooth random walk (see generate_target).
%   Shared across all three cases via p.seed_target for fair comparison.
%
% CONTROL (per active agent i):
%   u_i = v_T  +  k_r*K_i*(r*rhat_i - rho_i)  +  k_th*(errL+errR)*that_i
%
% RING: bidirectional, fixed at topology change, based on angular sort.
%
% COOPERATIVE RECOVERY (when UAV3 rejoins):
%   C1 — info sharing: neighbors give UAV3 p_T and v_T
%   C2 — slot guidance: UAV3 drives to angular midpoint of nearest 2
%   C3 — gap opening: neighbors temporarily yield angular space
%
% ONLINE BETA (if p.online = true):
%   gradient = Var[E(t)]  →  beta(k+1) = clip(beta - eta*Var, -1, 1)

%% ── Unpack ───────────────────────────────────────────────────────────────
dt      = p.dt;
N_steps = p.N_steps;
t       = p.t;
n_UAV   = p.n_UAV;
online  = p.online;
eta     = 0; if online, eta = p.eta; end

r_tol  = 0.5;    % orbit proximity threshold: Phase A → B [m]
k_slot = 3.0;    % slot guidance gain (C2)
k_coop = 0.6;    % gap opening gain (C3)

%% ── Pre-generate target trajectory ──────────────────────────────────────
% Use fixed seed so all cases follow the same random path.
rng(p.seed_target);
[tgt_pos, tgt_vel] = generate_target(p);

%% ── Storage ──────────────────────────────────────────────────────────────
X_hist      = zeros(N_steps, n_UAV, 2);
Theta_hist  = zeros(N_steps, n_UAV);
Target_hist = zeros(N_steps, 2);
Eerr_hist   = zeros(N_steps, 1);
Spread_hist = zeros(N_steps, 1);
K_hist      = zeros(N_steps, n_UAV);
Sigma_hist  = zeros(N_steps, n_UAV);
Active_hist = zeros(N_steps, n_UAV);
Phase_hist  = zeros(N_steps, 1);
Beta_hist   = zeros(N_steps, 1);
Var_hist    = zeros(N_steps, 1);

%% ── Process noise ────────────────────────────────────────────────────────
rng(p.seed_noise);
base_noise = randn(N_steps, n_UAV, 2);

%% ── Initial conditions ───────────────────────────────────────────────────
pT = tgt_pos(1,:);
rng(p.seed_ic);
theta_init = (0:n_UAV-1)' * (2*pi/n_UAV);
x = zeros(n_UAV, 2);
for i = 1:n_UAV
    x(i,:) = pT + p.r0 * [cos(theta_init(i)), sin(theta_init(i))];
end
x     = x + 0.4*(rand(n_UAV,2) - 0.5);
theta = theta_init;
Sigma = 0.3*ones(n_UAV,1);

%% ── Orbit radius ─────────────────────────────────────────────────────────
beta    = p.beta;
r_orbit = max(1.5, p.r0 * (1 - beta * 0.3));
if online, r_orbit = p.r0; end   % fixed for online to avoid moving reference

%% ── Online beta state ────────────────────────────────────────────────────
E_ema   = 0.5;
var_ema = 0.1;
gamma_e = 0.10;
gamma_v = 0.05;

%% ── State machine ────────────────────────────────────────────────────────
% 0=normal, 1=dead, 2=phaseA(off orbit), 3=phaseB(on orbit)
uav3_state = 0;
slot_pos   = zeros(1,2);
nb_L = 0; nb_R = 0;

%% ── Initial ring (4-UAV diamond) ─────────────────────────────────────────
active   = 1:n_UAV;
Delta    = 2*pi/4;
[LN, RN] = make_ring(active, theta);

fprintf('  Orbit r=%.2fm  beta_init=%+.2f  online=%d\n', r_orbit, beta, online);

%% ══════════════════════════════════════════════════════════════════════════
%% Simulation loop
%% ══════════════════════════════════════════════════════════════════════════
for k = 1:N_steps
    pT = tgt_pos(k,:);
    vT = tgt_vel(k,:);

    %% ── State transitions ────────────────────────────────────────────────
    if comm_loss

        % 0 → 1: UAV3 loses comms
        if uav3_state == 0 && t(k) >= p.t_fail
            uav3_state = 1;
            active     = [1, 2, 4];
            Delta      = 2*pi/3;
            [LN, RN]   = make_ring(active, theta);
            fprintf('  [t=%5.1fs]  UAV3 LOST → triangle\n', t(k));
        end

        % 1 → 2 or 3: rejoin beacon received
        if uav3_state == 1 && t(k) >= p.t_rejoin
            % C1: UAV3 receives pT from neighbors
            rel      = x(3,:) - pT;
            theta(3) = atan2(rel(2), rel(1));
            dist3    = norm(rel);

            % Neighbor discovery: Euclidean nearest two
            d_to_3 = arrayfun(@(j) norm(x(j,:) - x(3,:)), [1,2,4]);
            [~, si] = sort(d_to_3);
            pool = [1,2,4];
            nb_L = pool(si(1));
            nb_R = pool(si(2));
            fprintf('  [t=%5.1fs]  UAV3 BACK  neighbors: U%d U%d\n', t(k), nb_L, nb_R);

            % C2: slot = angular midpoint of nb_L and nb_R on orbit
            th_mid   = theta(nb_L) + 0.5*wrap_to_pi(theta(nb_R) - theta(nb_L));
            slot_pos = pT + r_orbit*[cos(th_mid), sin(th_mid)];

            % Restore 4-agent ring
            active   = 1:n_UAV;
            Delta    = 2*pi/4;
            [LN, RN] = make_ring(active, theta);
            Sigma(3) = norm(x(3,:) - slot_pos)*0.4 + 0.1;

            if abs(dist3 - r_orbit) < r_tol
                uav3_state = 3;
                fprintf('              On orbit → Phase B\n');
            else
                uav3_state = 2;
                fprintf('              Off orbit (%.2fm) → Phase A\n', abs(dist3-r_orbit));
            end
        end

        % 2 → 3: UAV3 reaches orbit
        if uav3_state == 2
            dist3 = norm(x(3,:) - pT);
            if abs(dist3 - r_orbit) < r_tol
                uav3_state = 3;
                fprintf('  [t=%5.1fs]  UAV3 on orbit → Phase B\n', t(k));
            else
                % Update slot to track moving target
                th_mid   = theta(nb_L) + 0.5*wrap_to_pi(theta(nb_R) - theta(nb_L));
                slot_pos = pT + r_orbit*[cos(th_mid), sin(th_mid)];
            end
        end
    end

    %% ── Online beta update ───────────────────────────────────────────────
    if online && k > 1
        E_prev  = Eerr_hist(k-1);
        E_ema   = (1-gamma_e)*E_ema   + gamma_e*E_prev;
        var_ema = (1-gamma_v)*var_ema + gamma_v*(E_prev - E_ema)^2;
        beta    = max(-1.0, min(1.0, beta - eta*var_ema));
    end

    %% ── Gains ────────────────────────────────────────────────────────────
    K = compute_beta_gains(Sigma, beta, p.k_base, p.k_min, p.k_max, n_UAV);

    %% ── C3: gap opening ──────────────────────────────────────────────────
    coop_u = zeros(n_UAV, 2);
    if comm_loss && (uav3_state == 2 || uav3_state == 3) && nb_L > 0
        th3       = theta(3);
        th_slot   = atan2(slot_pos(2)-pT(2), slot_pos(1)-pT(1));
        ang_dist3 = abs(wrap_to_pi(th3 - th_slot));
        fade      = min(1, ang_dist3/(pi/6));
        for j = [nb_L, nb_R]
            rho_j  = x(j,:) - pT;
            dist_j = max(norm(rho_j), 1e-6);
            rhat_j = rho_j/dist_j;
            that_j = [-rhat_j(2), rhat_j(1)];
            gap_dir = sign(wrap_to_pi(th_slot - theta(j)));
            coop_u(j,:) = k_coop * fade * gap_dir * that_j;
        end
    end

    %% ── Control inputs ───────────────────────────────────────────────────
    u = zeros(n_UAV, 2);

    for i = active
        rho_i  = x(i,:) - pT;
        dist_i = max(norm(rho_i), 1e-6);
        rhat_i = rho_i/dist_i;
        that_i = [-rhat_i(2), rhat_i(1)];

        u_radial = p.k_radial * K(i) * (r_orbit*rhat_i - rho_i);
        u_track  = p.k_track * vT;

        if i == 3 && uav3_state == 2
            % Phase A: radial + slot guidance, no angular consensus
            u_slot = k_slot * K(3) * (slot_pos - x(3,:));
            u(3,:) = u_track + u_radial + u_slot;

        elseif i == 3 && uav3_state == 3
            % Phase B: full bidirectional consensus
            err_L     = wrap_to_pi(theta(3) - theta(LN(3)) - Delta);
            err_R     = wrap_to_pi(theta(3) - theta(RN(3)) + Delta);
            u_angular = -p.k_theta * (err_L + err_R) * that_i;
            u(3,:)    = u_track + u_radial + u_angular;

        elseif i ~= 3 || uav3_state == 0
            % Normal: bidirectional consensus + cooperative gap opening
            err_L     = wrap_to_pi(theta(i) - theta(LN(i)) - Delta);
            err_R     = wrap_to_pi(theta(i) - theta(RN(i)) + Delta);
            u_angular = -p.k_theta * (err_L + err_R) * that_i;
            u(i,:)    = u_track + u_radial + u_angular + coop_u(i,:);
        end
        % uav3_state == 1 (dead): u(3,:) stays zero → hover
    end

    %% ── Integrate ────────────────────────────────────────────────────────
    w = zeros(n_UAV, 2);
    for i = 1:n_UAV
        w(i,:) = p.sigma_process * squeeze(base_noise(k,i,:))';
    end
    x = x + dt*(u + w);

    %% ── Update theta ─────────────────────────────────────────────────────
    for i = active
        rel      = x(i,:) - pT;
        theta(i) = atan2(rel(2), rel(1));
    end

    %% ── Update sigma ─────────────────────────────────────────────────────
    for i = active
        e_i      = abs(norm(x(i,:) - pT) - r_orbit);
        Sigma(i) = sqrt((1-p.ema_alpha)*Sigma(i)^2 + p.ema_alpha*e_i^2);
    end

    %% ── Metrics ──────────────────────────────────────────────────────────
    r_errors      = arrayfun(@(i) abs(norm(x(i,:)-pT)-r_orbit), active);
    Eerr_hist(k)  = mean(r_errors);

    th_s           = sort(mod(theta(active), 2*pi));
    gaps           = diff([th_s; th_s(1)+2*pi]);
    Spread_hist(k) = std(gaps);

    %% ── Store ────────────────────────────────────────────────────────────
    X_hist(k,:,:)     = x;
    Theta_hist(k,:)   = theta';
    Target_hist(k,:)  = pT;
    K_hist(k,:)       = K';
    Sigma_hist(k,:)   = Sigma';
    Active_hist(k,:)  = ismember(1:n_UAV, active);
    Phase_hist(k)     = uav3_state;
    Beta_hist(k)      = beta;
    Var_hist(k)       = var_ema;
end

fprintf('  Final beta = %+.4f\n', beta);

%% ── Pack results ─────────────────────────────────────────────────────────
results.X_hist      = X_hist;
results.Theta_hist  = Theta_hist;
results.Target_hist = Target_hist;
results.Eerr_hist   = Eerr_hist;
results.Spread_hist = Spread_hist;
results.K_hist      = K_hist;
results.Sigma_hist  = Sigma_hist;
results.Active_hist = Active_hist;
results.Phase_hist  = Phase_hist;
results.Beta_hist   = Beta_hist;
results.Var_hist    = Var_hist;
results.r_orbit     = r_orbit;
results.comm_loss   = comm_loss;
results.online      = online;
end

%% ════════════════════════════════════════════════════════════════════════
function [LN, RN] = make_ring(active, theta)
% Sort active agents by angle, assign bidirectional ring neighbors.
% Called once per topology change — NOT every step.
n  = 4;
LN = zeros(1,n); RN = zeros(1,n);
N  = length(active);
if N < 2, return; end
[~, si] = sort(mod(theta(active), 2*pi));
ring    = active(si);
for ii = 1:N
    i     = ring(ii);
    LN(i) = ring(mod(ii-2, N)+1);
    RN(i) = ring(mod(ii,   N)+1);
end
end

%% ════════════════════════════════════════════════════════════════════════
function [pos, vel] = generate_target(p)
% GENERATE_TARGET  Pre-compute target trajectory for all timesteps.
%
% Branches on p.target_mode:
%   'lissajous' — analytic figure-8 (uses target_Ax/Ay/wx/wy)
%   'random'    — smooth random walk (uses target_speed/tau/bounds/seed)
%
% Returns pos (N×2) and vel (N×2) arrays so the sim loop just indexes
% into them — no function calls inside the hot loop.

N  = p.N_steps;
dt = p.dt;
t  = p.t;

if strcmp(p.target_mode, 'lissajous')
    %% ── Lissajous figure-8 ───────────────────────────────────────────────
    pos = [p.target_Ax*sin(p.target_wx*t'), ...
           p.target_Ay*sin(p.target_wy*t')];
    vel = [p.target_Ax*p.target_wx*cos(p.target_wx*t'), ...
           p.target_Ay*p.target_wy*cos(p.target_wy*t')];

else
    %% ── Smooth random walk ───────────────────────────────────────────────
    % Algorithm:
    %   acc(k) = (1-alpha)*acc(k-1) + alpha*noise   [low-pass filtered]
    %   vel(k) = vel(k-1) + dt*(acc + restore)       [integrate, capped]
    %   pos(k) = pos(k-1) + dt*vel(k)                [integrate]
    %
    % Soft boundary: restoring force beyond target_bounds keeps target
    % in view without hard walls.
    spd = p.target_speed;
    tau = p.target_tau;
    bnd = p.target_bounds;

    alpha     = dt / tau;
    noise_std = spd * sqrt(2/tau);

    pos = zeros(N, 2);
    vel = zeros(N, 2);
    acc = [0, 0];

    for k = 2:N
        acc = (1-alpha)*acc + alpha*noise_std*randn(1,2);

        d = norm(pos(k-1,:));
        if d > bnd
            restore = -0.8*(d-bnd)/d * pos(k-1,:);
        else
            restore = [0, 0];
        end

        vel(k,:) = vel(k-1,:) + dt*(acc + restore);

        v_now = norm(vel(k,:));
        if v_now > spd
            vel(k,:) = vel(k,:) * (spd/v_now);
        end

        pos(k,:) = pos(k-1,:) + dt*vel(k,:);
    end
end
end

%% ════════════════════════════════════════════════════════════════════════
function a = wrap_to_pi(a)
    a = mod(a + pi, 2*pi) - pi;
end