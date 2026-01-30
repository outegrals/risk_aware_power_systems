A1 = [0.5 1; -3 -1];
A2 = [-1 -3; 1 0.5];

tspan = [0 100];

% Grid for vector fields
[x1g, x2g] = meshgrid(-3:0.5:3, -3:0.5:3);

% Vector field for A1: xdot = A1*x
dx1_A1 = A1(1,1)*x1g + A1(1,2)*x2g;
dx2_A1 = A1(2,1)*x1g + A1(2,2)*x2g;

% Vector field for A2: xdot = A2*x
dx1_A2 = A2(1,1)*x1g + A2(1,2)*x2g;
dx2_A2 = A2(2,1)*x1g + A2(2,2)*x2g;

figure
hold on
grid on
axis equal
xlabel('x_1')
ylabel('x_2')
title('Overlaid Phase Portraits: \dot{x}=A_1x and \dot{x}=A_2x')

% Plot both vector fields (different styles)
q1 = quiver(x1g, x2g, dx1_A1, dx2_A1, 0.7, 'k');   % A1 in black
q2 = quiver(x1g, x2g, dx1_A2, dx2_A2, 0.7, 'r');   % A2 in red
q2.LineStyle = '--';                               % dashed arrows for A2

% Multiple initial conditions (same for both)
ICs = [1 1; -1 0; -2 1];


for k = 1:size(ICs,1)
    scatter(ICs(k,1), ICs(k,2), 50, 'filled', 'DisplayName','Initial Condition')
    % Trajectory for A1
    [~, xA1] = ode45(@(t,x) A1*x, tspan, ICs(k,:)');
    plot(xA1(:,1), xA1(:,2), 'b', 'LineWidth', 1.5, 'DisplayName','A1 Trajectory')

    % Trajectory for A2
    [~, xA2] = ode45(@(t,x) A2*x, tspan, ICs(k,:)');
    plot(xA2(:,1), xA2(:,2), 'm', 'LineWidth', 1.5, 'DisplayName','A2 Trajectory')
end

legend
