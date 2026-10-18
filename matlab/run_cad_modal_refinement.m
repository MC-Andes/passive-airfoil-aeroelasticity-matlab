function report = run_cad_modal_refinement(tangentPath)
%RUN_CAD_MODAL_REFINEMENT Retained modes with identical physical fluid mesh.
s=load(tangentPath,'model','fluid');m=s.model;f=s.fluid;time=(0:f.wakeRows-2)*f.dt;Y=cell(4,1);traces=Y;rows=[];counts=unique([1,3,6,f.modes]);
saved=load(fullfile(fileparts(tangentPath),'dynamic.mat'),'dynamics');assert(numel(saved.dynamics.time)==numel(time) && max(abs(saved.dynamics.time-time))<1e-10,'cad:Time','Baseline response does not match this tangent.');
for k=1:numel(counts)
 [model,fluid]=fsi.truncate(m,f,counts(k));q0=zeros(counts(k),1);q0(1)=1e-5/max(abs(model.observation(:,1)));
 if counts(k)==f.modes && norm(saved.dynamics.Q(:,1)-q0)<1e-14 && norm(saved.dynamics.V(:,1))<1e-14,r=saved.dynamics;else,r=fsi.linear_dynamic(model,fluid,time,q0,zeros(size(q0)));end
 Y{k}=model.observation*r.Q;fit=cad.decay_fit(time,r.Q(1,:));traces{k}=struct('dynamics',r,'fit',fit);rows(end+1,:)=[counts(k),fit.frequencyHz,fit.growthPerSecond,r.relativeEnergyBalance]; %#ok<AGROW>
end
errors=zeros(size(counts));for k=1:numel(counts),errors(k)=max(abs(Y{k}-Y{end}),[],'all')/1e-5;end
report=array2table([rows,errors'],'VariableNames',{'modes','frequency_Hz','growth_per_s','energy_defect','maximum_sensor_error_over_initial_amplitude'});output=fullfile(fileparts(fileparts(tangentPath)),'modal_refinement');if ~isfolder(output),mkdir(output);end;writetable(report,fullfile(output,'comparison.csv'));save(fullfile(output,'responses.mat'),'Y','traces','time','report','-v7.3');disp(report);assert(errors(end-1)<.02,'cad:Modes','Retained-mode refinement exceeds two percent.');
end
