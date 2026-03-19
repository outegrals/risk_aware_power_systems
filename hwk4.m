%% EEL 6683 HW4 Problem 1(c)
%  Cooperative circular formation tracking under switching topologies
%  Three agents, double integrator, 120-degree spacing on circle
clear; clc; close all;

%% Parameters
r     = 5;          % circle radius
omega = 0.5;        % angular velocity (rad/s)
alpha = 2;          % position gain
beta  = 2;          % pinning gain
gamma = 2;          % velocity coupling ratio

% Desired phase offsets (120 degrees apart)
phi = [0, 2*pi/3, 4*pi/3];

% Simulation
dt   = 0.01;
T    = 40;
time = 0:dt:T;
N    = length(time);

%% Define switching topologies (each is {L, B})
% Topology 1: 1->2, 1 pinned
L1 = [0 0 0; -1 1 0; 0 0 0];
B1 = diag([1, 0, 0]);

% Topology 2: 2->3, 2 pinned
L2 = [0 0 0; 0 0 0; 0 -1 1];
B2 = diag([0, 1, 0]);

% Topology 3: 3->1, 3 pinned
L3 = [0 0 0; 0 0 0; -1 0 1];
B3 = diag([0, 0, 1]);

% Union of these three is connected with all agents reachable from leader
% Switch every T_switch seconds
T_switch = 2.0;  % switching period

%% Virtual leader references
q0_ref = @(t,ph) r * [cos(omega*t + ph); sin(omega*t + ph)];
v0_ref = @(t,ph) r*omega * [-sin(omega*t + ph); cos(omega*t + ph)];
a0_ref = @(t,ph) -r*omega^2 * [cos(omega*t + ph); sin(omega*t + ph)];

%% Initial conditions (random offsets)
rng(42);
q = 4*randn(2,3);      % positions [2 x 3]
v = 0.5*randn(2,3);    % velocities [2 x 3]

%% Storage
q_hist = zeros(2, 3, N);
v_hist = zeros(2, 3, N);
err_hist = zeros(N, 3);
topo_hist = zeros(N, 1);

%% Simulation loop
for k = 1:N
    t = time(k);
    
    % Select topology based on time (cycle through 3 topologies)
    topo_idx = mod(floor(t / T_switch), 3) + 1;
    topo_hist(k) = topo_idx;
    
    switch topo_idx
        case 1, L = L1; B = B1;
        case 2, L = L2; B = B2;
        case 3, L = L3; B = B3;
    end
    
    H = L + B;
    
    % Store states
    q_hist(:,:,k) = q;
    v_hist(:,:,k) = v;
    
    % Compute control for each agent
    u = zeros(2,3);
    for i = 1:3
        % Tracking errors
        q_tilde_i = q(:,i) - q0_ref(t, phi(i));
        v_tilde_i = v(:,i) - v0_ref(t, phi(i));
        
        % Consensus term
        consensus = zeros(2,1);
        for j = 1:3
            if i ~= j && L(i,j) < 0  % neighbor
                q_tilde_j = q(:,j) - q0_ref(t, phi(j));
                v_tilde_j = v(:,j) - v0_ref(t, phi(j));
                consensus = consensus + abs(L(i,j)) * ...
                    ((q_tilde_j - q_tilde_i) + gamma*(v_tilde_j - v_tilde_i));
            end
        end
        
        % Pinning term
        pinning = -B(i,i) * (q_tilde_i + gamma*v_tilde_i);
        
        % Feedforward + feedback
        u(:,i) = a0_ref(t, phi(i)) + alpha*consensus + beta*pinning;
        
        % Track position error norm
        err_hist(k,i) = norm(q_tilde_i);
    end
    
    % Euler integration
    if k < N
        q = q + v * dt;
        v = v + u * dt;
    end
end

%% Plotting
figure('Position', [100 100 1400 900]);

% --- Subplot 1: Trajectories in 2D ---
subplot(2,2,1); hold on; axis equal; grid on;
colors = lines(3);

% Draw desired circle
theta_c = linspace(0, 2*pi, 200);
plot(r*cos(theta_c), r*sin(theta_c), 'k--', 'LineWidth', 1);

for i = 1:3
    qi = squeeze(q_hist(:,i,:));
    plot(qi(1,:), qi(2,:), '-', 'Color', colors(i,:), 'LineWidth', 1.2);
    plot(qi(1,1), qi(2,1), 'o', 'Color', colors(i,:), 'MarkerSize', 8, 'MarkerFaceColor', colors(i,:));
    plot(qi(1,end), qi(2,end), 's', 'Color', colors(i,:), 'MarkerSize', 10, 'MarkerFaceColor', colors(i,:));
    
    % Desired final position
    qd = q0_ref(time(end), phi(i));
    plot(qd(1), qd(2), 'p', 'Color', colors(i,:), 'MarkerSize', 14, 'MarkerFaceColor', 'none', 'LineWidth', 2);
end
title('Agent Trajectories (o=start, square=end, star=desired)');
xlabel('x'); ylabel('y');
legend('Circle', 'Agent 1', '', '', '', 'Agent 2', '', '', '', 'Agent 3', 'Location', 'best');

% --- Subplot 2: Tracking errors ---
subplot(2,2,2); hold on; grid on;
for i = 1:3
    plot(time, err_hist(:,i), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
end
xlabel('Time (s)'); ylabel('||q_i - q_0^i||');
title('Position Tracking Errors');
legend('Agent 1', 'Agent 2', 'Agent 3');

% --- Subplot 3: Active topology ---
subplot(2,2,3); hold on; grid on;
stairs(time, topo_hist, 'k-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Topology Index');
title('Switching Topology (cycles through 3 disconnected graphs)');
yticks([1 2 3]);
yticklabels({'1\rightarrow2, pin 1', '2\rightarrow3, pin 2', '3\rightarrow1, pin 3'});
ylim([0.5 3.5]);

% --- Subplot 4: Final formation snapshot ---
subplot(2,2,4); hold on; axis equal; grid on;
t_snap = time(end);
for i = 1:3
    qi_end = q_hist(:,i,end);
    qd_end = q0_ref(t_snap, phi(i));
    plot(qi_end(1), qi_end(2), 'o', 'Color', colors(i,:), ...
        'MarkerSize', 12, 'MarkerFaceColor', colors(i,:));
    plot(qd_end(1), qd_end(2), 'p', 'Color', colors(i,:), ...
        'MarkerSize', 16, 'LineWidth', 2);
end
plot(r*cos(theta_c), r*sin(theta_c), 'k--', 'LineWidth', 1);
plot(0, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
title(sprintf('Final Formation at t = %.0f s', t_snap));
xlabel('x'); ylabel('y');
legend('Agent 1', 'Desired 1', 'Agent 2', 'Desired 2', 'Agent 3', 'Desired 3', 'Location', 'best');

sgtitle('HW4 Problem 1: Cooperative Circular Formation under Switching Topologies');

fprintf('Final tracking errors:\n');
for i = 1:3
    fprintf('  Agent %d: %.6f\n', i, err_hist(end,i));
end