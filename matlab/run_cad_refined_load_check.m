function report = run_cad_refined_load_check(bundlePath,staticStates,solutionFolder,output)
%RUN_CAD_REFINED_LOAD_CHECK Full FEM for rigid and converged deformed VLM loads.
root=fileparts(fileparts(mfilename('fullpath')));s=load(bundlePath,'bundle');b=s.bundle;m=b.model;c=b.config;
if nargin<4,output=fullfile(root,'results/cad_original/final_load_validation');end
out=output;if ~isfolder(out),mkdir(out);end;cad.provenance(root,c,out);
history=load(staticStates,'profiles','kind');r=history.profiles{end}.equilibrium;kind='horseshoe';if isfield(history,'kind'),kind=history.kind;end
if strcmp(kind,'uvlm'),options=struct('wakeLength',1000*b.flow.referenceChord,'coreRadius',1e-5);m.aero=@(X,V,f)fsi.uvlm_static(X,V,f,options);end
assert(strcmp(history.profiles{end}.source,bundlePath),'cad:Load','Static state and requested bundle differ.');
rigid=m.aero(m.X0,zeros(numel(m.X0)/3,3),b.flow);Fr=b.surface.H'*rigid.nodalForces;Fe=b.surface.H'*r.aero.nodalForces;F=[reshape(Fr',[],1),reshape(Fe',[],1)];
mesh=fem3d.elevate(fem3d.read_mesh(fullfile(root,c.discretization.meshes{end}),c.geometry.axisOrder,c.geometry.origin_m));system=fem3d.assemble(mesh,fem3d.material_reference(string(c.material.reference)));
old=load(fullfile(solutionFolder,'solution.mat'),'result');response=fem3d.constraint_response(system,F,old.result.bc,sparse(0,size(system.K,1)));Y=b.sensor.O*response.base;
one=m.K\(m.PhiSurface'*reshape(rigid.nodalForces',[],1));reduced=[m.observation*one,m.observation*r.q];errors=zeros(1,2);for k=1:2,errors(k)=norm(Y(:,k)-reduced(:,k),inf)/max(norm(Y(:,k),inf),realmin);end
report=struct('aerodynamicMesh',[b.surface.nx,b.surface.ny],'aerodynamicModel',kind,'fullReducedSensorRelative',errors,'maximumFullSensors_m',max(abs(Y),[],1), ...
 'energyResidual',response.energyResidual,'scope','two full FEM right-hand sides: refined rigid exterior VLM and load at its converged deformed equilibrium; original connected CAD and linear K; no viscous pressure claim');
cad.write_json(fullfile(out,'comparison.json'),report);save(fullfile(out,'solutions.mat'),'Y','reduced','response','report','-v7.3');disp(report);assert(max(errors)<.02,'cad:Reduction','Refined load requires enrichment of the solid basis.');
end
