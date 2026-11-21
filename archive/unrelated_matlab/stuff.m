clc; clear; close all;

%% Canvas Setup
fig = figure('Color', [0.05 0.10 0.25], 'Position', [100 100 1800 500]);
ax = axes('Parent', fig);
set(ax, 'Color', [0.05 0.10 0.25], 'XColor', 'none', 'YColor', 'none');
hold on;
axis off;
xlim([0 110]);
ylim([-4 4]);

lw = 4;

%% =========================================================
%% GRID LINES
%% =========================================================
for x = 0:2:110
    plot([x x], [-4 4], 'Color', [0.10 0.15 0.35], 'LineWidth', 0.5);
end
for y = -4:0.5:4
    plot([0 110], [y y], 'Color', [0.10 0.15 0.35], 'LineWidth', 0.5);
end
for x = 0:10:110
    plot([x x], [-4 4], 'Color', [0.15 0.22 0.45], 'LineWidth', 1.0);
end
for y = -4:2:4
    plot([0 110], [y y], 'Color', [0.15 0.22 0.45], 'LineWidth', 1.0);
end
plot([0 110], [0 0], 'Color', [0.20 0.30 0.60], 'LineWidth', 1.5);
plot([55 55], [-4 4], 'Color', [0.20 0.30 0.60], 'LineWidth', 1.5);

%% =========================================================
%% HELPER FUNCTION
%% =========================================================
function draw_stroke(x0, x1, y0, y1, freq, amp, func, col, lw)
    t  = linspace(0, 1, 500);
    xp = x0 + t*(x1 - x0);
    yb = y0 + t*(y1 - y0);
    switch func
        case 'sin'
            yw = amp * sin(freq * t * 2 * pi);
        case 'cos'
            yw = amp * cos(freq * t * 2 * pi);
        case 'sq'
            yw = amp * sign(sin(freq * t * 2 * pi));
        case 'tri'
            yw = amp * asin(sin(freq * t * 2 * pi)) / (pi/2);
        case 'poly2'
            yw = amp * (4*t.^2 - 4*t);
        case 'poly3'
            yw = amp * (8*t.^3 - 12*t.^2 + 4*t);
        case 'exp'
            yw = amp * exp(-4*t) .* sin(freq * t * 2 * pi);
        case 'flat'
            yw = zeros(size(t));
    end
    plot(xp, yb + yw, 'Color', col, 'LineWidth', lw);
end

%% =========================================================
%% COLOR PALETTE
%% =========================================================
gold   = [1.00 0.84 0.20];
red    = [0.90 0.30 0.30];
green  = [0.40 0.85 0.55];
blue   = [0.40 0.70 1.00];
pink   = [1.00 0.55 0.80];
orange = [0.95 0.55 0.20];
purple = [0.70 0.40 1.00];

%% =========================================================
%% LETTERS
%% =========================================================

%% === W ===
draw_stroke(2,  5,  2.5, -2.0, 1, 0.18, 'sin', gold, lw);
draw_stroke(5,  8,  -2.0, 0.5, 1, 0.18, 'sin', gold, lw);
draw_stroke(8,  11, 0.5, -2.0, 1, 0.18, 'sin', gold, lw);
draw_stroke(11, 14, -2.0, 2.5, 1, 0.18, 'sin', gold, lw);

%% === E ===
draw_stroke(16, 16, -2.0, 2.5,  3, 0.15, 'cos', red, lw);
draw_stroke(16, 22, 2.5,  2.5,  2, 0.15, 'sq',  red, lw);
draw_stroke(16, 21, 0.2,  0.2,  2, 0.15, 'sq',  red, lw);
draw_stroke(16, 22, -2.0, -2.0, 2, 0.15, 'sq',  red, lw);

%% === L ===
draw_stroke(24, 24, -2.0, 2.5,  3, 0.15, 'exp', green, lw);
draw_stroke(24, 30, -2.0, -2.0, 2, 0.15, 'tri', green, lw);

%% === C ===
t  = linspace(0.35*pi, 1.65*pi, 800);
xC = 36.5 + 3.2*cos(t);
yC = 0.25 + 3.2*sin(t);
ripple = 0.15 * sin(linspace(0, 6*pi, 800));
nx = -sin(t); ny = cos(t);
plot(xC + ripple.*nx, yC + ripple.*ny, 'Color', blue, 'LineWidth', lw);

%% === O ===
t  = linspace(0, 2*pi, 1000);
xO = 44.5 + 3.2*cos(t);
yO = 0.25 + 2.8*sin(t);
ripple = 0.15 * sin(linspace(0, 8*pi, 1000));
nx = -sin(t); ny = cos(t);
plot(xO + ripple.*nx, yO + ripple.*ny, 'Color', pink, 'LineWidth', lw);

%% === M ===
draw_stroke(48, 48, -2.0, 2.5,  3, 0.15, 'sin',   orange, lw);
draw_stroke(48, 51, 2.5,  0.0,  2, 0.15, 'poly2', orange, lw);
draw_stroke(51, 54, 0.0,  2.5,  2, 0.15, 'poly2', orange, lw);
draw_stroke(54, 54, -2.0, 2.5,  3, 0.15, 'sin',   orange, lw);

%% === E ===
draw_stroke(56, 56, -2.0, 2.5,  4, 0.12, 'sin', purple, lw);
draw_stroke(56, 62, 2.5,  2.5,  3, 0.12, 'tri', purple, lw);
draw_stroke(56, 60, 0.2,  0.2,  3, 0.12, 'tri', purple, lw);
draw_stroke(56, 62, -2.0, -2.0, 3, 0.12, 'tri', purple, lw);

%% =========================================================
%% SMILEY FACE
%% =========================================================
cx = 72;   % center x
cy = 0.25; % center y
r  = 3.0;  % radius
smiley_color = [1.00 0.84 0.20]; % gold

%% --- Outer circle (face outline) ---
t = linspace(0, 2*pi, 1000);
rx = 8;  % wider x radius
ry = 3.0;  % keep y radius same
ripple = 0.10 * sin(linspace(0, 10*pi, 1000));
nx = -sin(t); ny = cos(t);
xFace = cx + rx*cos(t);
yFace = cy + ry*sin(t);
plot(xFace + ripple.*nx, yFace + ripple.*ny, ...
    'Color', smiley_color, 'LineWidth', lw);

% --- Left eye (small upward parabola arc) ---
t_eye = linspace(0, pi, 200);
xEyeL = (cx - 1.1) + 0.6*cos(t_eye);
yEyeL = (cy + 1.0) + 0.35*sin(t_eye);
plot(xEyeL, yEyeL, 'Color', smiley_color, 'LineWidth', lw);

% --- Right eye ---
xEyeR = (cx + 1.1) + 0.6*cos(t_eye);
yEyeR = (cy + 1.0) + 0.35*sin(t_eye);
plot(xEyeR, yEyeR, 'Color', smiley_color, 'LineWidth', lw);

% --- Smile (downward arc using lower half ellipse) ---
t_smile = linspace(0.2*pi, 0.8*pi, 400);
xSmile  = cx + 1.6*cos(t_smile);
ySmile  = (cy - 0.3) - 1.0*sin(t_smile); % flipped to curve downward = smile
% add sine ripple on smile
ripple_s = 0.08 * sin(linspace(0, 4*pi, 400));
plot(xSmile, ySmile + ripple_s, 'Color', smiley_color, 'LineWidth', lw);

%% =========================================================
%% EXPORT
%% =========================================================
exportgraphics(fig, 'welcome_banner.png', ...
    'Resolution', 300, ...
    'BackgroundColor', [0.05 0.10 0.25]);

disp('Banner exported: welcome_banner.png');