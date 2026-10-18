function report = run_cad_static_convergence(bundlePaths,output,kind)
%RUN_CAD_STATIC_CONVERGENCE Same physical wind and solid, finer material lattices.
if nargin<2,output=fullfile(fileparts(fileparts(mfilename('fullpath'))),'results/cad_original/static_high_resolution');end
if nargin<3,kind='horseshoe';end;assert(ismember(string(kind),["horseshoe","uvlm"]),'cad:Options','Unknown steady aerodynamic model.');
if ~isfolder(output),mkdir(output);end;rows=[];previous=[];profiles=cell(numel(bundlePaths),1);
for k=1:numel(bundlePaths)
 s=load(bundlePaths{k},'bundle');b=s.bundle;m=b.model;f=b.flow;
 if strcmp(kind,'uvlm'),options=struct('wakeLength',1000*f.referenceChord,'coreRadius',1e-5);m.aero=@(X,V,flow)fsi.uvlm_static(X,V,flow,options);end
 if k==1,cad.provenance(fileparts(fileparts(mfilename('fullpath'))),b.config,output);end
 a=m.aero(m.X0,zeros(numel(m.X0)/3,3),f);one=m.K\(m.PhiSurface'*reshape(a.nodalForces',[],1));
 r=fsi.static(m,f,struct('initialQ',one,'relaxation',1,'loadTolerance',1e-7,'displacementTolerance',1e-7));Y=m.observation*r.q;
 change=NaN;if ~isempty(previous),change=norm(Y-previous,inf)/max(norm(Y,inf),realmin);end;previous=Y;
 rows(end+1,:)=[b.surface.nx,b.surface.ny,a.CL,r.aero.CL,max(abs(m.observation*one)),max(abs(Y)),change,r.relativeEquilibriumResidual,r.iterations]; %#ok<AGROW>
 profiles{k}=struct('source',bundlePaths{k},'eta',b.sensor.eta,'TE',Y,'equilibrium',r);
 fprintf('Static %dx%d: TE %.6g mm, change %.3g\n',b.surface.nx,b.surface.ny,1e3*max(abs(Y)),change);
 report=array2table(rows,'VariableNames',{'nx','ny','rigid_CL','coupled_CL','one_way_TE_m','two_way_TE_m','sensor_vector_relative_change','residual','iterations'});report.flow_model=repmat(string(kind),height(report),1);writetable(report,fullfile(output,'comparison.csv'));save(fullfile(output,'states.mat'),'profiles','report','kind','-v7.3');
end
end
