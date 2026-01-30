A1 = [0.5 1; -3 -1];
A2 = [-1 -3; 1 0.5];

tspan = [0 10];
Ts = 0.5;   % switching period (change as you like)

% Grid for vector fields (still useful for context)
[x1g, x2g] = meshgrid(-3:0.5:3, -3:0.5:3);

% Vector field for A1
dx1_A1 = A1(1,1)*x1g + A1(1,2)*x2g;
dx2_A1 = A1(2,1)*x1g + A1(2,2)*x2g;

% Vector field for A2
dx1_A2 = A2(1,1)*x1g + A2(1,2)*x2g;
dx2_A2 = A2(2,1)*x1g + A2(2,2)*x2g;

% Switched dynamics: alternate A1 and A2 every Ts
odefun_sw = @(t,x) ( mod(floor(t/Ts),2)==0 ) * (A1*x) + ...
                   ( mod(floor(t/Ts),2)==1 ) * (A2*x);

figure
hold on
grid on
axis equal
xlabel('x_1')
ylabel('x_2')

% Plot both vector fields (optional, but you had it)
q1 = quiver(x1g, x2g, dx1_A1, dx2_A1, 0.7, 'k');   % A1 in black
q2 = quiver(x1g, x2g, dx1_A2, dx2_A2, 0.7, 'r');   % A2 in red
q2.LineStyle = '--';

% Multiple initial conditions
ICs = [1 1; -1 0; -2 1];

% Plot switched trajectories
for k = 1:size(ICs,1)
    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [~, xS] = ode45(odefun_sw, tspan, ICs(k,:)', opts);
    plot(xS(:,1), xS(:,2), 'LineWidth', 1.8, 'DisplayName', 'Alternating Trajectories')

end

