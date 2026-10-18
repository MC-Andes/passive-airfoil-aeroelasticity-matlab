function report = run_cad_static_resolution(bundlePaths)
%RUN_CAD_STATIC_RESOLUTION Independent panel/wake and small-load static limits.
root=fileparts(fileparts(mfilename('fullpath')));rows=[];profiles={};
for k=1:numel(bundlePaths)
 s=load(bundlePaths{k},'bundle');b=s.bundle;m=b.model;f=b.flow;equilibrium=fsi.static(m,f,struct('relaxation',.6,'loadTolerance',1e-7,'displacementTolerance',1e-7));
 for length=[20,100,200]
  a=pae.vlm(m.X0,zeros(numel(m.X0)/3,3),f,struct('wakeLength',length*f.referenceChord));
  panels=pae.panels(m.X0);q=m.K\(m.PhiSurface'*reshape((panels.Hforce'*a.panelForces)',[],1));
  rows(end+1,:)=[b.surface.nx,b.surface.ny,length,a.CL,a.Cm,max(abs(m.observation*q)),equilibrium.aero.CL,max(abs(m.observation*equilibrium.q))]; %#ok<AGROW>
 end
 profiles{k}=struct('nx',b.surface.nx,'eta',b.sensor.eta,'TE',m.observation*equilibrium.q);
end
report=array2table(rows,'VariableNames',{'nx','ny','wake_chords','rigid_CL','rigid_Cm','one_way_TE_m','two_way_CL','two_way_TE_m'});output=fullfile(root,'results/cad_original/aerodynamic_resolution');if ~isfolder(output),mkdir(output);end
writetable(report,fullfile(output,'panel_wake.csv'));save(fullfile(output,'profiles.mat'),'profiles');
limits=[];
for U=[.01,.1,1,5,10]
 f.speed=U;a=m.aero(m.X0,zeros(numel(m.X0)/3,3),f);one=m.K\(m.PhiSurface'*reshape(a.nodalForces',[],1));two=fsi.static(m,f,struct('relaxation',.6,'loadTolerance',1e-8,'displacementTolerance',1e-8));limits(end+1,:)=[U,norm(m.observation*(two.q-one),inf)/max(norm(m.observation*one,inf),1e-20)]; %#ok<AGROW>
end
writematrix(limits,fullfile(output,'small_load_limit.csv'));assert(limits(1,2)<1e-5,'cad:FSI','One-way small-load limit failed.');
end
