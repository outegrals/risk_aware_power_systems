function plot_beta_results(t, dt, beta, X_hist, Eerr_hist, K_hist, ...
    Sigma_hist, Jtot_hist, formation_pos, n_UAV, t_fail, t_rejoin)
% PLOT_BETA_RESULTS  Two light-theme figures for a single beta run.
%
% Figure 1:  Trajectories (left)  +  Per-UAV x(t)/y(t) (right)
% Figure 2:  Formation error / All-UAV gains / All-UAV sigmas

%% ── Light-theme palette ──────────────────────────────────────────────────
BG       = [1.00 1.00 1.00];
AX_FG    = [0.15 0.15 0.15];
GRID_C   = [0.88 0.88 0.88];
DEAD_C   = [1.00 0.89 0.89];
FAIL_C   = [0.84 0.10 0.10];
REJOIN_C = [0.13 0.60 0.13];
REF_DASH = [0.65 0.65 0.65];

uav_colors = [0.13 0.47 0.71;   % UAV 1  steel blue
              0.20 0.63 0.17;   % UAV 2  green
              0.89 0.10 0.11;   % UAV 3  red
              0.99 0.55 0.10];  % UAV 4  orange
uav_names  = {'UAV 1','UAV 2','UAV 3','UAV 4'};

if     beta < -0.6,  col = [0.84 0.19 0.15];
elseif beta < -0.1,  col = [0.93 0.52 0.09];
elseif beta <=  0.1, col = [0.20 0.20 0.20];
elseif beta <=  0.6, col = [0.12 0.47 0.71];
else,                col = [0.09 0.25 0.60];
end

%% ── Key indices ──────────────────────────────────────────────────────────
[~, k_fail]   = min(abs(t - t_fail));
[~, k_rejoin] = min(abs(t - t_rejoin));
[~, k_ss]     = min(abs(t - (t(end) - 3.0)));

dead_px  = [t_fail t_rejoin t_rejoin t_fail];
beta_str = sprintf('\\beta = %.2f', beta);

    function style_ax(ax)
        set(ax,'Color',BG,'XColor',AX_FG,'YColor',AX_FG, ...
            'GridColor',GRID_C,'MinorGridColor',GRID_C,'GridAlpha',1.0, ...
            'FontSize',10,'FontName','Helvetica', ...
            'TickDir','out','Box','off','TitleFontWeight','bold');
    end

%% ════════════════════════════════════════════════════════════════════════
%  FIGURE 1 — Trajectories  +  Per-UAV state
% ════════════════════════════════════════════════════════════════════════
fig1 = figure('Name',sprintf('Formation — beta=%.2f',beta), ...
              'Color',BG,'Position',[40 100 1260 560]);

%% ── Left: 2D trajectories ────────────────────────────────────────────────
ax1 = subplot(1,2,1);
style_ax(ax1);
hold(ax1,'on'); grid(ax1,'on'); axis(ax1,'equal');
title(ax1,sprintf('UAV Trajectories  (%s)',beta_str),'Color',AX_FG);
xlabel(ax1,'X [m]','Color',AX_FG);
ylabel(ax1,'Y [m]','Color',AX_FG);

% Dead window shade — draw first so it's behind everything
patch(ax1, dead_px,[-6 -6 6 6],DEAD_C, ...
      'FaceAlpha',1.0,'EdgeColor','none','HandleVisibility','off');

% Desired diamond
des = [formation_pos; formation_pos(1,:)];
plot(ax1,des(:,1),des(:,2),'--','Color',REF_DASH, ...
     'LineWidth',1.2,'DisplayName','Desired slot');
scatter(ax1,formation_pos(:,1),formation_pos(:,2),90,REF_DASH,'x', ...
        'LineWidth',2.0,'HandleVisibility','off');
for i = 1:n_UAV
    text(ax1,formation_pos(i,1)+0.20,formation_pos(i,2)+0.20, ...
         uav_names{i},'FontSize',7.5,'Color',[0.50 0.50 0.50]);
end

% UAV trajectories
for i = 1:n_UAV
    tx = squeeze(X_hist(:,i,1));
    ty = squeeze(X_hist(:,i,2));
    uc = uav_colors(i,:);
    faded = uc + (BG - uc)*0.55;

    % Phase 1: nominal
    plot(ax1,tx(1:k_fail),ty(1:k_fail),'-', ...
         'Color',uc,'LineWidth',2.0,'HandleVisibility','off');

    % Phase 2: drift/degraded
    if i == 3
        plot(ax1,tx(k_fail:k_rejoin),ty(k_fail:k_rejoin),':', ...
             'Color',uc,'LineWidth',2.0,'HandleVisibility','off');
    else
        plot(ax1,tx(k_fail:k_rejoin),ty(k_fail:k_rejoin),'--', ...
             'Color',faded,'LineWidth',1.2,'HandleVisibility','off');
    end

    % Phase 3: recovery + SS
    plot(ax1,tx(k_rejoin:end),ty(k_rejoin:end),'-', ...
         'Color',uc,'LineWidth',2.0,'HandleVisibility','off');

    % Start marker — filled circle, in legend
    plot(ax1,tx(1),ty(1),'o', ...
         'Color',uc,'MarkerFaceColor',uc,'MarkerSize',8, ...
         'DisplayName',sprintf('%s  start',uav_names{i}));

    % End marker — open square, in legend
    plot(ax1,tx(end),ty(end),'s', ...
         'Color',uc,'MarkerFaceColor',BG,'MarkerEdgeColor',uc, ...
         'LineWidth',1.8,'MarkerSize',8, ...
         'DisplayName',sprintf('%s  end',uav_names{i}));
end

% UAV 3 events
plot(ax1,X_hist(k_fail,  3,1),X_hist(k_fail,  3,2),'x', ...
     'Color',FAIL_C,'MarkerSize',13,'LineWidth',2.5, ...
     'DisplayName','UAV3  comms lost');
plot(ax1,X_hist(k_rejoin,3,1),X_hist(k_rejoin,3,2),'^', ...
     'Color',REJOIN_C,'MarkerFaceColor',REJOIN_C,'MarkerSize',9, ...
     'DisplayName','UAV3  comms back');

legend(ax1,'Location','eastoutside','FontSize',7.5, ...
       'TextColor',AX_FG,'Color',BG,'EdgeColor',GRID_C);
xlim(ax1,[-6 6]); ylim(ax1,[-6 6]);

%% ── Right: Per-UAV x(t) and y(t) ────────────────────────────────────────
ax2 = subplot(1,2,2);
style_ax(ax2);
hold(ax2,'on'); grid(ax2,'on');
title(ax2,sprintf('Per-UAV State Over Time  (%s)',beta_str),'Color',AX_FG);
xlabel(ax2,'Time [s]','Color',AX_FG);
ylabel(ax2,'Position [m]','Color',AX_FG);

dim_style  = {'-','--'};
dim_labels = {'x','y'};

for dim = 1:2
    for i = 1:n_UAV
        vals = squeeze(X_hist(:,i,dim));
        plot(ax2,t,vals,dim_style{dim}, ...
             'Color',uav_colors(i,:),'LineWidth',1.6, ...
             'DisplayName',sprintf('%s (%s)',uav_names{i},dim_labels{dim}));

        yline(ax2,formation_pos(i,dim),':', ...
              'Color',uav_colors(i,:)+(BG-uav_colors(i,:))*0.60, ...
              'LineWidth',0.9,'HandleVisibility','off');
    end
end

xline(ax2,t_fail,  '--','Color',FAIL_C,  'LineWidth',1.2, ...
      'Label','Comms lost','HandleVisibility','off', ...
      'LabelOrientation','horizontal','LabelColor',FAIL_C);
xline(ax2,t_rejoin,'--','Color',REJOIN_C,'LineWidth',1.2, ...
      'Label','Comms back','HandleVisibility','off', ...
      'LabelOrientation','horizontal','LabelHorizontalAlignment','right', ...
      'LabelColor',REJOIN_C);

yl = ylim(ax2);
patch(ax2,dead_px,[yl(1) yl(1) yl(2) yl(2)],DEAD_C, ...
      'FaceAlpha',1.0,'EdgeColor','none','HandleVisibility','off');

legend(ax2,'Location','eastoutside','FontSize',7,'NumColumns',2, ...
       'TextColor',AX_FG,'Color',BG,'EdgeColor',GRID_C);
xlim(ax2,[t(1) t(end)]);

sgtitle(fig1,sprintf('Figure 1 — Trajectories & State    \\beta = %.2f',beta), ...
        'FontSize',13,'FontWeight','bold','Color',AX_FG);


%% ════════════════════════════════════════════════════════════════════════
%  FIGURE 2 — Error / All-UAV gains / All-UAV sigmas
% ════════════════════════════════════════════════════════════════════════
fig2 = figure('Name',sprintf('Signals — beta=%.2f',beta), ...
              'Color',BG,'Position',[80 60 960 780]);

%% ── Top: Formation error E(t) ────────────────────────────────────────────
ax3 = subplot(3,1,1);
style_ax(ax3);
hold(ax3,'on'); grid(ax3,'on');
title(ax3,'Mean Formation Error  E(t)','Color',AX_FG);
ylabel(ax3,'Error [m]','Color',AX_FG);

plot(ax3,t,Eerr_hist,'-','Color',col,'LineWidth',2.2);

ss_mean = mean(Eerr_hist(k_ss:end));
yline(ax3,ss_mean,'--','Color',col,'LineWidth',1.2, ...
      'Label',sprintf('SS = %.3f m',ss_mean), ...
      'LabelHorizontalAlignment','left','FontSize',8,'LabelColor',col);

xline(ax3,t_fail,  '--','Color',FAIL_C,  'LineWidth',1.2, ...
      'Label','Comms lost','HandleVisibility','off', ...
      'LabelOrientation','horizontal','LabelColor',FAIL_C);
xline(ax3,t_rejoin,'--','Color',REJOIN_C,'LineWidth',1.2, ...
      'Label','Comms back','HandleVisibility','off', ...
      'LabelOrientation','horizontal','LabelHorizontalAlignment','right', ...
      'LabelColor',REJOIN_C);

yl = ylim(ax3);
patch(ax3,dead_px,[yl(1) yl(1) yl(2) yl(2)],DEAD_C, ...
      'FaceAlpha',1.0,'EdgeColor','none','HandleVisibility','off');
xlim(ax3,[t(1) t(end)]);
set(ax3,'XTickLabel',[]);

%% ── Middle: All-UAV gains k_i(t) ────────────────────────────────────────
ax4 = subplot(3,1,2);
style_ax(ax4);
hold(ax4,'on'); grid(ax4,'on');
title(ax4, ...
    sprintf('Risk-Adjusted Gains  k_i(t)  [\\beta = %.2f drives these]',beta), ...
    'Color',AX_FG);
ylabel(ax4,'Gain  k_i','Color',AX_FG);

for i = 1:n_UAV
    ki = K_hist(:,i);
    % Blank UAV 3 during dead window — it has no active gain
    if i == 3
        ki(k_fail:k_rejoin) = NaN;
    end
    plot(ax4,t,ki,'-','Color',uav_colors(i,:),'LineWidth',1.8, ...
         'DisplayName',uav_names{i});
end

% Reference line at k_base
yline(ax4,1.5,':', 'Color',[0.55 0.55 0.55],'LineWidth',1.0, ...
      'Label','k_{base}','LabelHorizontalAlignment','left', ...
      'FontSize',8,'LabelColor',[0.55 0.55 0.55]);

xline(ax4,t_fail,  '--','Color',FAIL_C,  'LineWidth',1.2,'HandleVisibility','off');
xline(ax4,t_rejoin,'--','Color',REJOIN_C,'LineWidth',1.2,'HandleVisibility','off');

yl = ylim(ax4);
patch(ax4,dead_px,[yl(1) yl(1) yl(2) yl(2)],DEAD_C, ...
      'FaceAlpha',1.0,'EdgeColor','none','HandleVisibility','off');

legend(ax4,'Location','northeast','FontSize',9, ...
       'TextColor',AX_FG,'Color',BG,'EdgeColor',GRID_C);
xlim(ax4,[t(1) t(end)]);
set(ax4,'XTickLabel',[]);

%% ── Bottom: All-UAV uncertainty Sigma_i(t) ──────────────────────────────
ax5 = subplot(3,1,3);
style_ax(ax5);
hold(ax5,'on'); grid(ax5,'on');
title(ax5,'Tracking Uncertainty  \Sigma_i(t)  [input to gain law]', ...
      'Color',AX_FG);
xlabel(ax5,'Time [s]','Color',AX_FG);
ylabel(ax5,'\Sigma_i [m]','Color',AX_FG);

for i = 1:n_UAV
    si = Sigma_hist(:,i);
    if i == 3
        si(k_fail:k_rejoin) = NaN;
    end
    plot(ax5,t,si,'-','Color',uav_colors(i,:),'LineWidth',1.8, ...
         'DisplayName',uav_names{i});
end

xline(ax5,t_fail,  '--','Color',FAIL_C,  'LineWidth',1.2, ...
      'Label','Comms lost','HandleVisibility','off', ...
      'LabelOrientation','horizontal','LabelColor',FAIL_C);
xline(ax5,t_rejoin,'--','Color',REJOIN_C,'LineWidth',1.2, ...
      'Label','Comms back','HandleVisibility','off', ...
      'LabelOrientation','horizontal','LabelHorizontalAlignment','right', ...
      'LabelColor',REJOIN_C);

yl = ylim(ax5);
patch(ax5,dead_px,[yl(1) yl(1) yl(2) yl(2)],DEAD_C, ...
      'FaceAlpha',1.0,'EdgeColor','none','HandleVisibility','off');

legend(ax5,'Location','northeast','FontSize',9, ...
       'TextColor',AX_FG,'Color',BG,'EdgeColor',GRID_C);
xlim(ax5,[t(1) t(end)]);

sgtitle(fig2, ...
    sprintf('Figure 2 — Signals Over Time    \\beta = %.2f', beta), ...
    'FontSize',13,'FontWeight','bold','Color',AX_FG);

set(findall(fig1,'-property','FontName'),'FontName','Helvetica');
set(findall(fig2,'-property','FontName'),'FontName','Helvetica');

end