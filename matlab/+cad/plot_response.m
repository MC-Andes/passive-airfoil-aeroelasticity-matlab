function plot_response(sensor,pressures,Y,frf,config,output)
%PLOT_RESPONSE Simple export with explicit light theme and conditional labels.
fig=figure('Visible','off','Color','w','Position',[100,100,1200,480]);
if isprop(fig,'Theme'), fig.Theme='light'; end
tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
ax=nexttile; plot(ax,sensor.eta,1e3*Y,'-o','LineWidth',1); grid(ax,'on');
xlabel(ax,'\eta'); ylabel(ax,'w_z cerca de TE [mm]');
legend(ax,compose('%g Pa',pressures),'Location','best');
title(ax,'CAD: presión analítica impuesta; escenario A','FontSize',11);
ax=nexttile; semilogy(ax,frf.frequencyHz,max(abs(frf.Y),[],1),'LineWidth',1); grid(ax,'on');
xlabel(ax,'f [Hz]'); ylabel(ax,'máximo |w| [m/Pa]');
title(ax,sprintf('FRF modal: amortiguamiento supuesto %.1f%%',100*config.discretization.dampingRatio),'FontSize',11);
exportgraphics(fig,fullfile(output,'response.png'),'Resolution',150); close(fig);
end
