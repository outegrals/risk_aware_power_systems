function play_encirclement(p, r1, r2)
% PLAY_ENCIRCLEMENT  Interactive trajectory viewer with time slider.
%
% USAGE (after running main_encirclement):
%   play_encirclement(p, results1, results2)
%
% Two side-by-side panels — Case 1 (left) and Case 2 (right).
% Drag the slider to scrub through time. Press Play to animate.
%
% What you see at each timestep:
%   - Faded trail  : full path up to current time
%   - Solid dot    : current UAV position
%   - Orbit circle : desired orbit at current target position
%   - Target dot   : current target position
%   - Time needle  : vertical line on a mini error plot below

%% ── Palette ──────────────────────────────────────────────────────────────
BG       = [1.00 1.00 1.00];
AX_FG    = [0.15 0.15 0.15];
GRID_C   = [0.88 0.88 0.88];
DEAD_C   = [1.00 0.89 0.89];
RECOV_C  = [0.88 0.97 0.88];
FAIL_C   = [0.84 0.10 0.10];
REJOIN_C = [0.13 0.60 0.13];
TARGET_C = [0.50 0.00 0.80];

uav_colors = [0.13 0.47 0.71;
              0.20 0.63 0.17;
              0.89 0.10 0.11;
              0.99 0.55 0.10];
uav_names = {'UAV 1','UAV 2','UAV 3','UAV 4'};

t       = p.t;
N_steps = p.N_steps;
n_UAV   = p.n_UAV;
beta    = p.beta;

[~, k_fail]   = min(abs(t - p.t_fail));
[~, k_rejoin] = min(abs(t - p.t_rejoin));

%% ── Figure layout ────────────────────────────────────────────────────────
fig = figure('Name', sprintf('Formation Player  beta=%.2f', beta), ...
             'Color', BG, 'Position', [30 60 1400 860]);

% Two trajectory axes (top row)
ax_c1   = subplot('Position', [0.03  0.28  0.44  0.65]);  % Case 1 traj
ax_c2   = subplot('Position', [0.52  0.28  0.44  0.65]);  % Case 2 traj

% Two mini error axes (bottom row — show time context)
ax_e1   = subplot('Position', [0.03  0.10  0.44  0.13]);  % Case 1 error
ax_e2   = subplot('Position', [0.52  0.10  0.44  0.13]);  % Case 2 error

    function style_ax(ax)
        set(ax,'Color',BG,'XColor',AX_FG,'YColor',AX_FG, ...
            'GridColor',GRID_C,'GridAlpha',1.0, ...
            'FontSize',10,'FontName','Helvetica', ...
            'TickDir','out','Box','off','TitleFontWeight','bold');
    end

style_ax(ax_c1); style_ax(ax_c2);
style_ax(ax_e1); style_ax(ax_e2);

%% ── Static elements: error plots with full time series ───────────────────
% Draw once — only the time needle redraws
function draw_error_background(ax, Eerr, is_loss, label)
    hold(ax,'on'); grid(ax,'on');
    title(ax, label, 'Color',AX_FG,'FontSize',9);
    xlabel(ax,'Time [s]','Color',AX_FG);
    ylabel(ax,'Err [m]','Color',AX_FG,'FontSize',8);

    if is_loss
        yl_est = [0 max(Eerr)*1.2+0.01];
        patch(ax,[p.t_fail p.t_rejoin p.t_rejoin p.t_fail], ...
              [yl_est(1) yl_est(1) yl_est(2) yl_est(2)], DEAD_C, ...
              'FaceAlpha',1.0,'EdgeColor','none','HandleVisibility','off');
        patch(ax,[p.t_rejoin t(end) t(end) p.t_rejoin], ...
              [yl_est(1) yl_est(1) yl_est(2) yl_est(2)], RECOV_C, ...
              'FaceAlpha',0.6,'EdgeColor','none','HandleVisibility','off');
    end

    plot(ax, t, Eerr, '-','Color',[0.3 0.3 0.3],'LineWidth',1.4);
    xlim(ax,[t(1) t(end)]);
end

draw_error_background(ax_e1, r1.Eerr_hist, r1.comm_loss, 'Case 1  radial error');
draw_error_background(ax_e2, r2.Eerr_hist, r2.comm_loss, 'Case 2  radial error');

%% ── Slider and controls ──────────────────────────────────────────────────
% Time slider
uicontrol('Style','text','Units','normalized', ...
    'Position',[0.10 0.045 0.08 0.03], ...
    'String','Time [s]:','BackgroundColor',BG, ...
    'ForegroundColor',AX_FG,'FontSize',10,'HorizontalAlignment','right');

slider = uicontrol('Style','slider','Units','normalized', ...
    'Position',[0.19 0.048 0.55 0.025], ...
    'Min',1,'Max',N_steps,'Value',1, ...
    'SliderStep',[1/(N_steps-1) 10/(N_steps-1)]);

time_label = uicontrol('Style','text','Units','normalized', ...
    'Position',[0.75 0.045 0.08 0.03], ...
    'String','t = 0.00s','BackgroundColor',BG, ...
    'ForegroundColor',AX_FG,'FontSize',10,'HorizontalAlignment','left');

% Play / Pause button
is_playing = false;
play_btn = uicontrol('Style','pushbutton','Units','normalized', ...
    'Position',[0.84 0.042 0.07 0.038], ...
    'String','▶  Play','FontSize',10,'BackgroundColor',[0.94 0.97 1.0], ...
    'Callback', @toggle_play);

% Speed control
uicontrol('Style','text','Units','normalized', ...
    'Position',[0.92 0.045 0.04 0.03], ...
    'String','Speed:','BackgroundColor',BG, ...
    'ForegroundColor',AX_FG,'FontSize',9,'HorizontalAlignment','right');
speed_slider = uicontrol('Style','slider','Units','normalized', ...
    'Position',[0.93 0.048 0.065 0.022], ...
    'Min',1,'Max',20,'Value',5, ...
    'SliderStep',[1/19 5/19]);

%% ── Initial draw ─────────────────────────────────────────────────────────
draw_frame(1);
slider.Callback = @(src,~) draw_frame(round(src.Value));

%% ── Play/pause logic ─────────────────────────────────────────────────────
    function toggle_play(~,~)
        is_playing = ~is_playing;
        if is_playing
            play_btn.String = '⏸  Pause';
            play_btn.BackgroundColor = [1.0 0.95 0.90];
            animate();
        else
            play_btn.String = '▶  Play';
            play_btn.BackgroundColor = [0.94 0.97 1.0];
        end
    end

    function animate()
        k = round(slider.Value);
        step = max(1, round(speed_slider.Value));
        while is_playing && ishandle(fig)
            k = k + step;
            if k > N_steps
                k = 1;   % loop
            end
            slider.Value = k;
            draw_frame(k);
            drawnow limitrate;
        end
    end

%% ── Core draw function ───────────────────────────────────────────────────
    function draw_frame(k)
        tk = t(k);
        set(time_label,'String',sprintf('t = %.2fs', tk));

        draw_case_frame(ax_c1, ax_e1, r1, k, tk, 'Case 1  (no comms loss)');
        draw_case_frame(ax_c2, ax_e2, r2, k, tk, 'Case 2  (UAV3 comms loss)');
    end

    function draw_case_frame(ax_traj, ax_err, res, k, tk, label)
        X       = res.X_hist;
        Tgt     = res.Target_hist;
        r_orb   = res.r_orbit;
        is_loss = res.comm_loss;

        % Determine phase from Phase_hist
        phase_k = 0;
        if is_loss && isfield(res,'Phase_hist')
            phase_k = res.Phase_hist(k);
        end
        uav3_active = (phase_k ~= 1);

        % Phase colors
        PHASE_A_C = [0.98 0.95 0.80];   % yellow: guided approach
        PHASE_B_C = [0.88 0.97 0.88];   % green:  consensus recovery

        %% ── Trajectory axes ──────────────────────────────────────────────
        cla(ax_traj);
        hold(ax_traj,'on'); grid(ax_traj,'on'); axis(ax_traj,'equal');

        title(ax_traj, sprintf('%s    t = %.2fs', label, tk), ...
              'Color',AX_FG,'FontSize',11);
        xlabel(ax_traj,'X [m]','Color',AX_FG);
        ylabel(ax_traj,'Y [m]','Color',AX_FG);

        % Phase background
        if is_loss
            if phase_k == 1
                fill(ax_traj,[-9 9 9 -9],[-9 -9 9 9],DEAD_C, ...
                     'EdgeColor','none','FaceAlpha',0.5);
                text(ax_traj,-8,8,'UAV3 hovering  (triangle)', ...
                     'FontSize',9,'Color',FAIL_C,'FontWeight','bold');
            elseif phase_k == 2
                fill(ax_traj,[-9 9 9 -9],[-9 -9 9 9],PHASE_A_C, ...
                     'EdgeColor','none','FaceAlpha',0.7);
                text(ax_traj,-8,8,{'Cooperative recovery', 'Phase A: guided to orbit'}, ...
                     'FontSize',9,'Color',[0.70 0.50 0.00],'FontWeight','bold');
            elseif phase_k == 3
                fill(ax_traj,[-9 9 9 -9],[-9 -9 9 9],PHASE_B_C, ...
                     'EdgeColor','none','FaceAlpha',0.5);
                text(ax_traj,-8,8,{'Cooperative recovery', 'Phase B: angular consensus'}, ...
                     'FontSize',9,'Color',[0.10 0.55 0.10],'FontWeight','bold');
            end
        end

        % Target: full faded path + current position
        plot(ax_traj, Tgt(1:k,1), Tgt(1:k,2), '-', ...
             'Color',[TARGET_C 0.25],'LineWidth',1.2,'HandleVisibility','off');
        plot(ax_traj, Tgt(k,1), Tgt(k,2), 'p', ...
             'Color',TARGET_C,'MarkerFaceColor',TARGET_C,'MarkerSize',16, ...
             'DisplayName','Target');

        % Current orbit circle
        th = linspace(0,2*pi,200);
        plot(ax_traj, Tgt(k,1)+r_orb*cos(th), Tgt(k,2)+r_orb*sin(th), '--', ...
             'Color',[0.55 0.55 0.55 0.40],'LineWidth',1.0, ...
             'HandleVisibility','off');

        % UAV trails and current positions
        for i = 1:n_UAV
            uc = uav_colors(i,:);

            % Trail up to current time
            trail_x = squeeze(X(1:k,i,1));
            trail_y = squeeze(X(1:k,i,2));

            % For UAV 3 in loss case during dead window — dotted trail
            if is_loss && i == 3 && ~uav3_active
                plot(ax_traj, trail_x, trail_y, ':', ...
                     'Color',[uc 0.35],'LineWidth',1.5,'HandleVisibility','off');
            else
                plot(ax_traj, trail_x, trail_y, '-', ...
                     'Color',[uc 0.30],'LineWidth',1.5,'HandleVisibility','off');
            end

            % Current position — large filled circle
            if is_loss && i == 3 && ~uav3_active
                % UAV 3 hovering — show with dashed ring to indicate hold
                plot(ax_traj, X(k,i,1), X(k,i,2), 'o', ...
                     'Color',uc,'MarkerFaceColor',BG, ...
                     'MarkerEdgeColor',uc,'LineWidth',2.0,'MarkerSize',14, ...
                     'DisplayName',sprintf('%s (hover)',uav_names{i}));
            else
                plot(ax_traj, X(k,i,1), X(k,i,2), 'o', ...
                     'Color',uc,'MarkerFaceColor',uc,'MarkerSize',14, ...
                     'DisplayName',uav_names{i});
            end

            % UAV label next to current position
            text(ax_traj, X(k,i,1)+0.25, X(k,i,2)+0.25, ...
                 sprintf('U%d',i),'FontSize',8,'Color',uc,'FontWeight','bold');
        end

        % Draw formation lines between active UAVs
        if uav3_active
            active = 1:n_UAV;
        else
            active = [1,2,4];
        end
        n_act = length(active);
        for ii = 1:n_act
            i = active(ii);
            j = active(mod(ii,n_act)+1);
            plot(ax_traj,[X(k,i,1) X(k,j,1)],[X(k,i,2) X(k,j,2)],'-', ...
                 'Color',[0.70 0.70 0.70],'LineWidth',0.8, ...
                 'HandleVisibility','off');
        end

        legend(ax_traj,'Location','southeast','FontSize',8, ...
               'TextColor',AX_FG,'Color',BG,'EdgeColor',GRID_C);
        xlim(ax_traj,[-9 9]); ylim(ax_traj,[-9 9]);

        %% ── Time needle on error plot ─────────────────────────────────────
        % Delete old needle if it exists, draw new one
        delete(findobj(ax_err,'Tag','needle'));
        yl = ylim(ax_err);
        plot(ax_err,[tk tk],yl,'-','Color',[0.20 0.20 0.80], ...
             'LineWidth',1.8,'Tag','needle','HandleVisibility','off');
    end

end