function report = run_cad_load_validation(bundlePath,viscousFolder,structuralCampaign,basisFile)
%RUN_CAD_LOAD_VALIDATION Full own FEM vs reduction for independently refined loads.
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<3,structuralCampaign=fullfile(root,'results/cad_original/campaign_20261017T190836945');end
if nargin<4,basisFile=fullfile(root,'results/cad_original/aeroelastic_primary/basis.mat');end
b=load(bundlePath,'bundle');b=b.bundle;c=b.config;mesh=fem3d.elevate(fem3d.read_mesh(fullfile(root,string(c.discretization.meshes{end})),c.geometry.axisOrder,c.geometry.origin_m));system=fem3d.assemble(mesh,fem3d.material_reference(string(c.material.reference)));
aero=fsi.vlm(b.model.X0,zeros(numel(b.model.X0)/3,3),b.flow);Fa=b.surface.H'*aero.nodalForces;vc=jsondecode(fileread(fullfile(viscousFolder,'case.json')));rows=readmatrix(fullfile(viscousFolder,'surface_traction.csv'));viscous=cad.section_traction_transfer(mesh,c,rows,vc.Re,vc.alpha);Fv=viscous.nodalForces;
bc=struct('dofs',fem3d.dofs(find(mesh.nodes(:,1)<=c.supports.x_over_c*c.geometry.chord_m)));F=[full(reshape(Fa',[],1)),full(reshape(Fv',[],1))];response=fem3d.constraint_response(system,F,bc,sparse(0,size(system.K,1)));p=load(basisFile,'basis');sensor=cad.sensors(mesh,c);q=b.model.K\(p.basis'*F);fullSensors=sensor.O*response.base;reducedSensors=b.model.observation*q;relativeError=max(abs(fullSensors-reducedSensors),[],1)./max(abs(fullSensors),[],1);
guide=load(fullfile(structuralCampaign,'mechanisms/guide_compliance.mat'),'response');
for scenario=['B','C']
 selected=1:6;if scenario=='C',selected=2:6;end
 coordinates=guide.response.observation(selected,:);lambda=guide.response.compliance(selected,selected)\(coordinates*response.base(:,1));uGuide=response.base(:,1)-guide.response.influence(:,selected)*lambda;
 modeFile=fullfile(structuralCampaign,'mechanisms',"modes_"+scenario+"_level_3.mat");if ~isfile(modeFile),modeFile=fullfile(structuralCampaign,'mechanisms',"modes_"+scenario+'.mat');end
 modes=load(modeFile,'modes');
 cad.guided_bundle(system,modes.modes,uGuide,F(:,1),b,scenario,fullfile(root,'results/cad_original',"aeroelastic_guide_"+scenario));
end
report=struct('aerodynamicMesh',[b.surface.nx,b.surface.ny],'fullReducedSensorRelative',relativeError,'maximumFullSensors_m',max(abs(fullSensors),[],1),'energyResidual',response.energyResidual,'viscousSpeed_ms',viscous.speed_ms,'viscousForceDefect',viscous.forceDefect,'viscousMomentDefect',viscous.momentDefect,'scope','refined rigid VLM load and computed low-Re viscous section traction, separate full FEM solves; no turbulent pressure claim');output=fullfile(root,'results/cad_original/load_validation');if ~isfolder(output),mkdir(output);end;cad.write_json(fullfile(output,'comparison.json'),report);save(fullfile(output,'solutions.mat'),'fullSensors','reducedSensors','report','q','-v7.3');disp(report);assert(relativeError(1)<.02,'cad:Reduction','Refined aerodynamic load requires more structural modes.');
end
