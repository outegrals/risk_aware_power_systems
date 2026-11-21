function plot_encirclement(p, r1, r2, r3)
% PLOT_ENCIRCLEMENT
%   Figure 1 — Case 1 trajectory
%   Figure 2 — Case 2 trajectory
%   Figure 3 — Case 3 trajectory  (online beta)
%   Figure 4 — Case 2 signals
%   Figure 5 — Case 3 signals + beta(t) panel
%   Figure 6 — Comparison: Cases 2 vs 3 (fixed vs online beta)

%% ── Palette ──────────────────────────────────────────────────────────────
BG       = [1.00 1.00 1.00];
AX_FG    = [0.15 0.15 0.15];
GRID_C   = [0.88 0.88 0.88];
DEAD_C   = [1.00 0.89 0.89];
PHASE_AC = [0.98 0.95 0.80];
PHASE_BC = [0.88 0.97 0.88];
FAIL_C   = [0.84 0.10 0.10];
REJOIN_C = [0.13 0.60 0.13];
TARGET_C = [0.50 0.00 0.80];
C2_COL   = [0.84 0.19 0.15];   % Case 2 line color
C3_COL   = [0.12 0.47 0.71];   % Case 3 line color

uav_colors = [0.13 0.47 0.71;
              0.20 0.63 0.17;
              0.89 0.10 0.11;
              0.99 0.55 0.10];
uav_names = {'UAV 1','UAV 2','UAV 3','UAV 4'};

t     = p.t;
n_UAV = p.n_UAV;

[~, k_fail]   = min(abs(t - p.t_fail));
[~, k_rejoin] = min(abs(t - p.t_rejoin));
dead_px = [p.t_fail p.t_rejoin p.t_rejoin p.t_fail];

    function style_ax(ax)
        set(ax,'Color',BG,'XColor',AX_FG,'YColor',AX_FG, ...
            'GridColor',GRID_C,'GridAlpha',1.0, ...
            'FontSize',11,'FontName','Helvetica', ...
            'TickDir','out','Box','off','TitleFontWeight','bold');
    end

    function shade_phases(ax, res)
        if ~res.comm_loss, return; end
        yl = ylim(ax); if yl(1)==yl(2), yl=[0 1]; end
        yv = [yl(1) yl(1) yl(2) yl(2)];
        patch(ax,dead_px,yv,DEAD_C,'FaceAlpha',1.0,'EdgeColor','none','HandleVisibility','off');
        ph   = res.Phase_hist;
        idxA = find(ph==2); idxB = find(ph==3);
        if ~isempty(idxA)
            patch(ax,[t(idxA(1)) t(idxA(end)) t(idxA(end)) t(idxA(1))],yv, ...
                  PHASE_AC,'FaceAlpha',0.9,'EdgeColor','none','HandleVisibility','off');
        end
        if ~isempty(idxB)
            patch(ax,[t(idxB(1)) t(idxB(end)) t(idxB(end)) t(idxB(1))],yv, ...
                  PHASE_BC,'FaceAlpha',0.8,'EdgeColor','none','HandleVisibility','off');
        end
        xline(ax,p.t_fail,  '--','Color',FAIL_C,  'LineWidth',1.2,'Label','Lost', ...
              'HandleVisibility','off','LabelOrientation','horizontal','LabelColor',FAIL_C);
        xline(ax,p.t_rejoin,'--','Color',REJOIN_C,'LineWidth',1.2,'Label','Back', ...
              'HandleVisibility','off','LabelOrientation','horizontal', ...
              'LabelHorizontalAlignment','right','LabelColor',REJOIN_C);
    end

%% ── Trajectory helper ────────────────────────────────────────────────────
    function draw_traj(ax, res, ttl)
        X       = res.X_hist;
        Tgt     = res.Target_hist;
        r_orb   = res.r_orbit;
        is_loss = res.comm_loss;

        cla(ax); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
        title(ax,ttl,'Color',AX_FG,'FontSize',12);
        xlabel(ax,'X [m]','Color',AX_FG); ylabel(ax,'Y [m]','Color',AX_FG);

        plot(ax,Tgt(:,1),Tgt(:,2),'-','Color',[TARGET_C 0.30],'LineWidth',1.5,'DisplayName','Target path');
        plot(ax,Tgt(1,1),Tgt(1,2),'p','Color',TARGET_C,'MarkerFaceColor',TARGET_C,'MarkerSize',14,'DisplayName','Target start');
        plot(ax,Tgt(end,1),Tgt(end,2),'h','Color',TARGET_C,'MarkerFaceColor',TARGET_C,'MarkerSize',14,'DisplayName','Target end');

        th = linspace(0,2*pi,200);
        plot(ax,Tgt(1,1)+r_orb*cos(th),Tgt(1,2)+r_orb*sin(th),'--','Color',[0.6 0.6 0.6 0.25],'LineWidth',1.0,'HandleVisibility','off');
        plot(ax,Tgt(end,1)+r_orb*cos(th),Tgt(end,2)+r_orb*sin(th),'--','Color',[0.6 0.6 0.6 0.25],'LineWidth',1.0,'HandleVisibility','off');

        for i = 1:n_UAV
            tx = squeeze(X(:,i,1)); ty = squeeze(X(:,i,2));
            uc = uav_colors(i,:); fd = uc+(BG-uc)*0.55;
            if is_loss && i==3
                plot(ax,tx(1:k_fail),ty(1:k_fail),'-','Color',uc,'LineWidth',2.0,'HandleVisibility','off');
                plot(ax,tx(k_fail:k_rejoin),ty(k_fail:k_rejoin),'-.','Color',[uc 0.40],'LineWidth',1.4,'HandleVisibility','off');
                plot(ax,tx(k_rejoin:end),ty(k_rejoin:end),'-','Color',uc,'LineWidth',2.0,'HandleVisibility','off');
            elseif is_loss
                plot(ax,tx(1:k_fail),ty(1:k_fail),'-','Color',uc,'LineWidth',2.0,'HandleVisibility','off');
                plot(ax,tx(k_fail:k_rejoin),ty(k_fail:k_rejoin),'-','Color',fd,'LineWidth',1.3,'HandleVisibility','off');
                plot(ax,tx(k_rejoin:end),ty(k_rejoin:end),'-','Color',uc,'LineWidth',2.0,'HandleVisibility','off');
            else
                plot(ax,tx,ty,'-','Color',uc,'LineWidth',2.0,'HandleVisibility','off');
            end
            plot(ax,tx(1),ty(1),'o','Color',uc,'MarkerFaceColor',uc,'MarkerSize',9,'DisplayName',sprintf('%s start',uav_names{i}));
            plot(ax,tx(end),ty(end),'s','Color',uc,'MarkerFaceColor',BG,'MarkerEdgeColor',uc,'LineWidth',1.8,'MarkerSize',9,'DisplayName',sprintf('%s end',uav_names{i}));
        end

        if is_loss
            plot(ax,X(k_fail,3,1),X(k_fail,3,2),'x','Color',FAIL_C,'MarkerSize',15,'LineWidth',2.8,'DisplayName','UAV3 lost');
            plot(ax,X(k_rejoin,3,1),X(k_rejoin,3,2),'^','Color',REJOIN_C,'MarkerFaceColor',REJOIN_C,'MarkerSize',11,'DisplayName','UAV3 back');
        end

        legend(ax,'Location','eastoutside','FontSize',9,'TextColor',AX_FG,'Color',BG,'EdgeColor',GRID_C);
        xlim(ax,[-9 9]); ylim(ax,[-9 9]);
    end

%% ── Signal panels helper ─────────────────────────────────────────────────
    function draw_signals(fig, res, show_beta)
        if show_beta
            ax1 = subplot(4,1,1,'Parent',fig);
            ax2 = subplot(4,1,2,'Parent',fig);
            ax3 = subplot(4,1,3,'Parent',fig);
            ax4 = subplot(4,1,4,'Parent',fig);
        else
            ax1 = subplot(3,1,1,'Parent',fig);
            ax2 = subplot(3,1,2,'Parent',fig);
            ax3 = subplot(3,1,3,'Parent',fig);
            ax4 = [];
        end
        style_ax(ax1); style_ax(ax2); style_ax(ax3);

        % Radial error
        hold(ax1,'on'); grid(ax1,'on');
        title(ax1,'Mean Radial Error  E(t)','Color',AX_FG);
        ylabel(ax1,'Error [m]','Color',AX_FG);
        plot(ax1,t,res.Eerr_hist,'-','Color',[0.15 0.15 0.15],'LineWidth',2.0);
        shade_phases(ax1,res);
        xlim(ax1,[t(1) t(end)]); set(ax1,'XTickLabel',[]);

        % Angular spread
        hold(ax2,'on'); grid(ax2,'on');
        title(ax2,'Angular Spread  (std of gaps,  0 = perfect)','Color',AX_FG);
        ylabel(ax2,'Std [rad]','Color',AX_FG);
        plot(ax2,t,res.Spread_hist,'-','Color',[0.15 0.15 0.15],'LineWidth',2.0);
        shade_phases(ax2,res);
        xlim(ax2,[t(1) t(end)]); set(ax2,'XTickLabel',[]);

        % Gains
        hold(ax3,'on'); grid(ax3,'on');
        title(ax3,'Risk-Adjusted Gains  k_i(t)','Color',AX_FG);
        ylabel(ax3,'Gain  k_i','Color',AX_FG);
        for i = 1:n_UAV
            ki = res.K_hist(:,i);
            if res.comm_loss && i==3, ki(k_fail:k_rejoin)=NaN; end
            plot(ax3,t,ki,'-','Color',uav_colors(i,:),'LineWidth',1.8,'DisplayName',uav_names{i});
        end
        yline(ax3,p.k_base,':','Color',[0.55 0.55 0.55],'LineWidth',1.0,'Label','k_{base}','LabelHorizontalAlignment','left','FontSize',9,'LabelColor',[0.55 0.55 0.55]);
        shade_phases(ax3,res);
        legend(ax3,'Location','northeast','FontSize',9,'TextColor',AX_FG,'Color',BG,'EdgeColor',GRID_C);
        xlim(ax3,[t(1) t(end)]);

        % Beta(t) — online only
        if show_beta && ~isempty(ax4)
            style_ax(ax4);
            hold(ax4,'on'); grid(ax4,'on');
            title(ax4,'\beta(t) — online gradient descent','Color',AX_FG);
            xlabel(ax4,'Time [s]','Color',AX_FG);
            ylabel(ax4,'\beta','Color',AX_FG);
            plot(ax4,t,res.Beta_hist,'-','Color',C3_COL,'LineWidth',2.2);
            yline(ax4, 0,'--','Color',[0.5 0.5 0.5],'LineWidth',1.0,'Label','neutral','LabelColor',[0.5 0.5 0.5],'LabelHorizontalAlignment','left','FontSize',8);
            yline(ax4,-1,':','Color',FAIL_C,  'LineWidth',0.8,'HandleVisibility','off');
            yline(ax4, 1,':','Color',REJOIN_C,'LineWidth',0.8,'HandleVisibility','off');
            shade_phases(ax4,res);
            xlim(ax4,[t(1) t(end)]);
            ylim(ax4,[-1.1 1.1]);
            set(ax3,'XTickLabel',[]);
        end
    end

%% ════════════════════════════════════════════════════════════════════════
%  FIGURE 1 — Case 1 trajectory
% ════════════════════════════════════════════════════════════════════════
fig1 = figure('Name','Case 1 Trajectory','Color',BG,'Position',[20 120 780 720]);
ax = axes('Parent',fig1); style_ax(ax);
draw_traj(ax, r1, sprintf('Case 1 — No loss  (\\beta=%.2f)', p.beta));
set(findall(fig1,'-property','FontName'),'FontName','Helvetica');

%% ════════════════════════════════════════════════════════════════════════
%  FIGURE 2 — Case 2 trajectory (fixed beta)
% ════════════════════════════════════════════════════════════════════════
fig2 = figure('Name','Case 2 Trajectory','Color',BG,'Position',[820 120 780 720]);
ax = axes('Parent',fig2); style_ax(ax);
draw_traj(ax, r2, sprintf('Case 2 — Comms loss, fixed \\beta=%.2f', p.beta));
annotation(fig2,'textbox',[0.01 0.01 0.98 0.04],'String', ...
    'Pink=hovering   Yellow=Phase A (guided)   Green=Phase B (consensus)', ...
    'EdgeColor','none','FontSize',8,'Color',[0.4 0.4 0.4],'FitBoxToText','off');
set(findall(fig2,'-property','FontName'),'FontName','Helvetica');

%% ════════════════════════════════════════════════════════════════════════
%  FIGURE 3 — Case 3 trajectory (online beta)
% ════════════════════════════════════════════════════════════════════════
fig3 = figure('Name','Case 3 Trajectory','Color',BG,'Position',[20 80 780 720]);
ax = axes('Parent',fig3); style_ax(ax);
draw_traj(ax, r3, sprintf('Case 3 — Comms loss, online \\beta (init=%.2f, \\eta=%.3f)', ...
          r3.Beta_hist(1), p.eta));
annotation(fig3,'textbox',[0.01 0.01 0.98 0.04],'String', ...
    'Pink=hovering   Yellow=Phase A (guided)   Green=Phase B (consensus)', ...
    'EdgeColor','none','FontSize',8,'Color',[0.4 0.4 0.4],'FitBoxToText','off');
set(findall(fig3,'-property','FontName'),'FontName','Helvetica');

%% ════════════════════════════════════════════════════════════════════════
%  FIGURE 4 — Case 2 signals (fixed beta)
% ════════════════════════════════════════════════════════════════════════
fig4 = figure('Name','Case 2 Signals','Color',BG,'Position',[820 80 880 680]);
draw_signals(fig4, r2, false);
sgtitle(fig4,sprintf('Case 2 Signals — Fixed \\beta=%.2f',p.beta), ...
        'FontSize',13,'FontWeight','bold','Color',AX_FG);
set(findall(fig4,'-property','FontName'),'FontName','Helvetica');

%% ════════════════════════════════════════════════════════════════════════
%  FIGURE 5 — Case 3 signals + beta(t) (online beta)
% ════════════════════════════════════════════════════════════════════════
fig5 = figure('Name','Case 3 Signals + Beta','Color',BG,'Position',[60 40 880 800]);
draw_signals(fig5, r3, true);
sgtitle(fig5,sprintf('Case 3 Signals — Online \\beta  (\\eta=%.3f)',p.eta), ...
        'FontSize',13,'FontWeight','bold','Color',AX_FG);
set(findall(fig5,'-property','FontName'),'FontName','Helvetica');

%% ════════════════════════════════════════════════════════════════════════
%  FIGURE 6 — Comparison: Case 2 vs Case 3
% ════════════════════════════════════════════════════════════════════════
fig6 = figure('Name','Comparison: Fixed vs Online Beta','Color',BG,'Position',[100 30 1000 700]);

ax6a = subplot(3,1,1,'Parent',fig6); style_ax(ax6a);
hold(ax6a,'on'); grid(ax6a,'on');
title(ax6a,'Radial Error  E(t)','Color',AX_FG,'FontSize',11);
ylabel(ax6a,'Mean error [m]','Color',AX_FG);
plot(ax6a,t,r2.Eerr_hist,'-', 'Color',C2_COL,'LineWidth',2.0,'DisplayName',sprintf('Fixed \\beta=%.2f',p.beta));
plot(ax6a,t,r3.Eerr_hist,'--','Color',C3_COL,'LineWidth',2.0,'DisplayName','Online \beta');
shade_phases(ax6a,r2);
legend(ax6a,'Location','northeast','FontSize',10,'TextColor',AX_FG,'Color',BG,'EdgeColor',GRID_C);
xlim(ax6a,[t(1) t(end)]); set(ax6a,'XTickLabel',[]);

ax6b = subplot(3,1,2,'Parent',fig6); style_ax(ax6b);
hold(ax6b,'on'); grid(ax6b,'on');
title(ax6b,'Angular Spread  (lower = more uniform)','Color',AX_FG,'FontSize',11);
ylabel(ax6b,'Std of gaps [rad]','Color',AX_FG);
plot(ax6b,t,r2.Spread_hist,'-', 'Color',C2_COL,'LineWidth',2.0,'DisplayName',sprintf('Fixed \\beta=%.2f',p.beta));
plot(ax6b,t,r3.Spread_hist,'--','Color',C3_COL,'LineWidth',2.0,'DisplayName','Online \beta');
shade_phases(ax6b,r2);
legend(ax6b,'Location','northeast','FontSize',10,'TextColor',AX_FG,'Color',BG,'EdgeColor',GRID_C);
xlim(ax6b,[t(1) t(end)]); set(ax6b,'XTickLabel',[]);

ax6c = subplot(3,1,3,'Parent',fig6); style_ax(ax6c);
hold(ax6c,'on'); grid(ax6c,'on');
title(ax6c,'\beta(t) — online gradient descent','Color',AX_FG,'FontSize',11);
xlabel(ax6c,'Time [s]','Color',AX_FG);
ylabel(ax6c,'\beta','Color',AX_FG);
yline(ax6c, p.beta,'--','Color',C2_COL,'LineWidth',1.5,'Label',sprintf('Fixed \\beta=%.2f',p.beta),'LabelColor',C2_COL,'LabelHorizontalAlignment','left','FontSize',9);
plot(ax6c,t,r3.Beta_hist,'-','Color',C3_COL,'LineWidth',2.2,'DisplayName','Online \beta');
yline(ax6c, 0,'--','Color',[0.5 0.5 0.5],'LineWidth',0.8,'HandleVisibility','off');
shade_phases(ax6c,r2);
legend(ax6c,'Location','northeast','FontSize',10,'TextColor',AX_FG,'Color',BG,'EdgeColor',GRID_C);
xlim(ax6c,[t(1) t(end)]);
ylim(ax6c,[-1.1 1.1]);

% Annotate key events on beta plot
xline(ax6c,p.t_fail,  '--','Color',FAIL_C,  'LineWidth',1.2,'HandleVisibility','off','Label','Lost','LabelOrientation','horizontal','LabelColor',FAIL_C);
xline(ax6c,p.t_rejoin,'--','Color',REJOIN_C,'LineWidth',1.2,'HandleVisibility','off','Label','Back','LabelOrientation','horizontal','LabelHorizontalAlignment','right','LabelColor',REJOIN_C);

sgtitle(fig6,sprintf('Case 2 (Fixed \\beta=%.2f)  vs  Case 3 (Online \\beta)',p.beta), ...
        'FontSize',13,'FontWeight','bold','Color',AX_FG);
set(findall(fig6,'-property','FontName'),'FontName','Helvetica');

end