function render_microgrid_evidence(out)
T=readtable(fullfile(out,'representative_trajectory.csv'));
f=figure('Visible','off','Color','w','Position',[0 0 1200 450]);tiledlayout(1,2);
nexttile;plot(T.hour_end,[T.min_soc T.mean_soc T.max_soc]);ylim([0 1]);xlabel('Hours');ylabel('State of charge');legend('Min','Mean','Max');grid on;
nexttile;plot(T.hour_end,[T.mean_forecast_MW T.mean_cap_MW T.RN_mean_bid_MW T.retuned_mean_bid_MW]);xlabel('Hours');ylabel('Mean power (MW)');legend('Forecast','Cap','RN','Retuned');grid on;
set(findall(f,'Type','axes'),'Color','w','XColor','k','YColor','k');set(findall(f,'Type','text'),'Color','k');set(findall(f,'Type','legend'),'Color','w','TextColor','k');
exportgraphics(f,fullfile(out,'microgrid_trajectory.pdf'),'ContentType','vector');exportgraphics(f,fullfile(out,'microgrid_trajectory.png'));close(f);
% Native TikZ export for manuscript inclusion; no source document is edited.
file=fopen(fullfile(out,'trajectory_plot.tex'),'w');clean=onCleanup(@()fclose(file));
fprintf(file,'\\begin{figure*}[t]\n\\centering\n');
keys={{'min_soc','mean_soc','max_soc'},{'mean_forecast_MW','mean_cap_MW','RN_mean_bid_MW','retuned_mean_bid_MW'}};colors={{'gray,dashed','blue,thick','gray,dashed'},{'black,thick','gray,dotted','orange,dashed','blue,thick'}};
for panel=1:2
 fprintf(file,'\\begin{tikzpicture}[x=0.145cm,y=%.5fcm,font=\\scriptsize]\n',3/panel);
 fprintf(file,'\\draw[->] (0,0)--(49,0) node[right]{h};\n\\draw[->] (0,0)--(0,%g);\n',panel);
 for tick=[0 12 24 36 48],fprintf(file,'\\node[below] at (%g,0){%g};\n',tick,tick);end
 for j=1:numel(keys{panel})
  fprintf(file,'\\draw[%s] plot coordinates {',colors{panel}{j});v=T.(keys{panel}{j});fprintf(file,'(%.2f,%.5f) ',[T.hour_end v]');fprintf(file,'};\n');
 end
 fprintf(file,'\\end{tikzpicture}\n');if panel==1,fprintf(file,'\\hspace{0.7cm}\n');end
end
fprintf(file,'\\caption{Sequential synthetic microgrid operation: battery state of charge (left) and forecasts, caps, and bids (right). Shared bid-independent physical policy.}\n\\end{figure*}\n');
end
