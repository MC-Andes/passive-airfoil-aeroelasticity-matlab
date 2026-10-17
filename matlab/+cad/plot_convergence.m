function plot_convergence(folder4,folder10,output)
%PLOT_CONVERGENCE Actual structural observables and residuals; no validity claim.
a=jsondecode(fileread(fullfile(folder4,'acceptance.json'))); b=jsondecode(fileread(fullfile(folder10,'acceptance.json')));
fig=figure('Visible','off','Color','w','Position',[100,100,1250,1000]);
if isprop(fig,'Theme'), fig.Theme='light'; end
tiledlayout(fig,3,2,'Padding','compact','TileSpacing','compact');
sets={a,b}; folders={folder4,folder10}; labels={'Tet4','Tet10'}; styles={'-o','-s'};
for j=1:2
    ax=nexttile; hold(ax,'on'); rows=sets{j}.summary;
    for k=1:numel(rows)
        data=readtable(fullfile(folders{j},rows(k).name,'sensors.csv'));
        plot(ax,data.eta,1e5*data.w_m_per_Pa,styles{j},'LineWidth',1.1, ...
            'DisplayName',sprintf('%d DOFs',3*rows(k).nodes));
    end
    grid(ax,'on'); xlabel(ax,'\eta'); ylabel(ax,'w_z cerca de TE [mm]');
    title(ax,[labels{j},': 100 Pa impuestos']); legend(ax,'Location','best');
end
ax=nexttile; hold(ax,'on');
for j=1:2
    rows=sets{j}.summary; loglog(ax,3*[rows.nodes],1e5*[rows.maxSensor_m_per_Pa],styles{j},'LineWidth',1.2,'DisplayName',labels{j});
end
ax.XScale='log'; ax.YScale='log'; grid(ax,'on'); xlabel(ax,'DOFs'); ylabel(ax,'máximo |w_z| a 100 Pa [mm]'); legend(ax,'Location','best');
ax=nexttile; hold(ax,'on');
for j=1:2
    rows=sets{j}.summary; freq=zeros(numel(rows),1);
    for k=1:numel(rows), modes=readtable(fullfile(folders{j},rows(k).name,'modes.csv')); freq(k)=modes.frequency_Hz(1); end
    semilogx(ax,3*[rows.nodes],freq,styles{j},'LineWidth',1.2,'DisplayName',labels{j});
end
ax.XScale='log'; grid(ax,'on'); xlabel(ax,'DOFs'); ylabel(ax,'primera frecuencia [Hz]'); legend(ax,'Location','best');
ax=nexttile; hold(ax,'on');
for j=1:2
    rows=sets{j}.summary; loglog(ax,3*[rows.nodes],[rows.residual],styles{j},'LineWidth',1.2,'DisplayName',labels{j});
end
ax.XScale='log'; ax.YScale='log'; yline(ax,1e-6,'--','Criterio inicial','DisplayName','criterio 10^{-6}'); grid(ax,'on'); xlabel(ax,'DOFs'); ylabel(ax,'residuo relativo original'); legend(ax,'Location','best');
ax=nexttile; hold(ax,'on');
for j=1:2
    rows=sets{j}.summary; loglog(ax,3*[rows.nodes],1e4*[rows.energy_J_per_Pa2],styles{j},'LineWidth',1.2,'DisplayName',labels{j});
end
ax.XScale='log'; ax.YScale='log'; grid(ax,'on'); xlabel(ax,'DOFs'); ylabel(ax,'energía a 100 Pa [J]'); legend(ax,'Location','best');
sgtitle(fig,'CAD original, escenario A: diagnóstico sin convergencia ni validación experimental');
exportgraphics(fig,output,'Resolution',160); close(fig);
end
