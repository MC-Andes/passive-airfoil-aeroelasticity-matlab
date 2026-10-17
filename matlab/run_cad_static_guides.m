function report = run_cad_static_guides()
%RUN_CAD_STATIC_GUIDES Actual-CAD A/B/C equilibrium map under the same wind.
% Flow, skin mapping and original connectivity are held fixed. The guide is
% a declared hypothesis and all pressures are exterior thin-surface loads.
root=fileparts(fileparts(mfilename('fullpath')));paths={'aeroelastic_128x26','aeroelastic_guide_B','aeroelastic_guide_C'};rows=[];
for scenario=1:3
 s=load(fullfile(root,'results/cad_original',paths{scenario},'reduced_cad.mat'),'bundle');b=s.bundle;m=b.model;f=b.flow;q0=zeros(size(m.K,1),1);
 for speed=[2.5,5,7.5,10]
  f.speed=speed;a=m.aero(m.X0,zeros(numel(m.X0)/3,3),f);one=m.K\(m.PhiSurface'*reshape(a.nodalForces',[],1));r=fsi.static(m,f,struct('initialQ',q0,'relaxation',1,'loadTolerance',1e-7,'displacementTolerance',1e-7));q0=r.q;
  rows(end+1,:)=[scenario,speed,a.CL,r.aero.CL,max(abs(m.observation*one)),max(abs(m.observation*r.q)),r.relativeEquilibriumResidual,r.iterations]; %#ok<AGROW>
 end
 slower=fsi.static(m,f,struct('relaxation',.6,'loadTolerance',1e-7,'displacementTolerance',1e-7));assert(norm(m.observation*(slower.q-r.q),inf)/max(norm(m.observation*r.q,inf),realmin)<1e-4,'cad:Guide','Guide equilibrium depends on relaxation.');
end
report=array2table(rows,'VariableNames',{'scenario_ABC','U_ms','rigid_CL','coupled_CL','one_way_TE_m','coupled_TE_m','residual','iterations'});out=fullfile(root,'results/cad_original/guide_flow');if ~isfolder(out),mkdir(out);end;writetable(report,fullfile(out,'static_map.csv'));cad.write_json(fullfile(out,'scope.json'),struct('cavityPressure','ambient through original open ports; no enclosed-air dynamics in this static map','scenarios','A monolithic, B six bar coordinates blocked, C chord translation allowed','experimentalValidation',false));disp(report);
fig=figure('Visible','off','Color','w');if isprop(fig,'Theme'),fig.Theme='light';end;hold on;for scenario=1:3,take=report.scenario_ABC==scenario;plot(report.U_ms(take),1e3*report.coupled_TE_m(take),'o-','DisplayName',char(64+scenario));end;xlabel('U [m/s]');ylabel('máximo |w_TE| [mm]');grid on;legend;exportgraphics(fig,fullfile(out,'guide_flow.png'),'Resolution',150);close(fig);
end
