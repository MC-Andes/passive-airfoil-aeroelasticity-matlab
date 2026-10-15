function fig = plot_study(model,summary,solutions,response)
%PLOT_STUDY Deformed mean surface, sweep and illustrative transient.
if nargin<4, response=[]; end
scale=1; yLabel='q [unidades indicadas en la leyenda]';
labels=string(model.structure.labels(:))+" ["+string(model.structure.coordinateUnits(:))+"]";
if all(strcmp(model.structure.coordinateUnits,'m'))
    scale=1e3; yLabel='q [mm]'; labels=string(model.structure.labels(:));
end
fig=figure('Color','w','Position',[80,80,1200,780]);
if isprop(fig,'Theme'), fig.Theme='light'; end
tiledlayout(2,2,'TileSpacing','compact');
nexttile; base=model.surface.nodes; deformed=solutions{end}.state.nodes;
surf(base(:,:,1),base(:,:,2),base(:,:,3),'FaceColor','none','EdgeColor',[.6,.6,.6]); hold on;
surf(deformed(:,:,1),deformed(:,:,2),deformed(:,:,3),'FaceColor',[.1,.5,.7],'FaceAlpha',.8);
axis equal; grid on; view(35,25); xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Superficie media: referencia y equilibrio');
nexttile; angles=unique(summary.alpha_deg,'stable'); hold on;
for angle=angles'
    rows=summary.alpha_deg==angle & logical(summary.converged);
    plot(summary.U_m_s(rows),summary.CL(rows),'-o','DisplayName',sprintf('%.1f grados',angle));
end
grid on; xlabel('U [m/s]'); ylabel('C_L'); legend('Location','best'); title('Barrido estático');
nexttile; hold on;
rows=summary.alpha_deg==angles(end) & logical(summary.converged);
for j=1:model.ndof
    name=summary.Properties.VariableNames{8+j};
    plot(summary.U_m_s(rows),scale*summary.(name)(rows),'-o','DisplayName',labels(j));
end
grid on; xlabel('U [m/s]'); ylabel(yLabel);
legend('Location','best','Interpreter','none'); title('Coordenadas estructurales');
nexttile;
if ~isempty(response)
    plot(response.time,scale*response.q'); grid on; xlabel('t [s]'); ylabel(yLabel);
legend(labels,'Location','best','Interpreter','none'); title('Transitorio cuasiestacionario');
else
    plot(summary.U_m_s,summary.residual,'o'); set(gca,'YScale','log'); grid on;
    xlabel('U [m/s]'); ylabel('Residual de equilibrio'); title('Convergencia');
end
sgtitle(string(model.structure.source),'Interpreter','none');
axesHandles=findall(fig,'Type','axes');
set(axesHandles,'Color','w','XColor',[.15,.15,.15],'YColor',[.15,.15,.15], ...
    'ZColor',[.15,.15,.15],'GridColor',[.65,.65,.65],'FontSize',10);
end
