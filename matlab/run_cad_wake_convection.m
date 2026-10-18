function report = run_cad_wake_convection(bundlePath)
%RUN_CAD_WAKE_CONVECTION Short rigid-CAD test of computed wake advection.
% This isolates convection; it is not a replacement for refined small-
% perturbation frozen-wake FSI or a claim about viscous vortex cores.
root=fileparts(fileparts(mfilename('fullpath')));s=load(bundlePath,'bundle');b=s.bundle;X=b.model.X0;velocity=zeros(numel(X)/3,3);flow=b.flow;cases=[0,1,.0005;1,1,.0005;1,2,.0005;1,2,.00025];rows=[];histories=cell(4,1);finalStates=histories;
output=fullfile(root,'results/cad_original/wake_convection');if ~isfolder(output),mkdir(output);end
for k=1:4
 convection='frozen';if cases(k,1),convection='free';end
 [state,base]=uvlm.initialize(X,velocity,flow,struct('convection',convection,'convectionOrder',cases(k,2),'wakeLength',1000*flow.referenceChord));time=0:cases(k,3):.05;trace=zeros(numel(time),4);trace(1,:)=[0,base.CL,base.Cm,0];
 for j=2:numel(time),[state,r]=uvlm.step(state,X,velocity,flow,cases(k,3));trace(j,:)=[time(j),r.CL,r.Cm,r.kelvinDefect];end
 straight=state.geometry.trailingEdge(:,3)'+(state.wake.nodes(2:end,:,1)-state.geometry.trailingEdge(:,1)')*tan(flow.alpha);deviation=state.wake.nodes(2:end,:,3)-straight;
 use=time>=.025;rows(end+1,:)=[k,cases(k,:),mean(trace(use,2)),trace(end,2),max(trace(:,4)),max(abs(deviation),[],'all')]; %#ok<AGROW>
 histories{k}=trace;finalStates{k}=state;
end
report=array2table(rows,'VariableNames',{'case','free_convection','order','dt_s','mean_CL_last_half','end_CL','kelvin_defect','maximum_wake_departure_from_ambient_line_m'});writetable(report,fullfile(output,'comparison.csv'));save(fullfile(output,'states.mat'),'histories','finalStates','report','-v7.3');cad.write_json(fullfile(output,'scope.json'),struct('rigidOriginalCADSurface',true,'duration_s',.05,'physicalModel','computed inviscid wake, frozen vs induced convection','limitations','short convection/control study on 8x6 CAD surface; refined main FSI uses its separately declared frozen-wake model; no viscous deficit inferred'));
fig=figure('Visible','off','Color','w');if isprop(fig,'Theme'),fig.Theme='light';end;hold on;for k=1:4,plot(histories{k}(:,1),histories{k}(:,2),'DisplayName',sprintf('caso %d',k));end;xlabel('t [s]');ylabel('CL calculado');legend;grid on;exportgraphics(fig,fullfile(output,'convection_loads.png'),'Resolution',150);close(fig);disp(report);
end
