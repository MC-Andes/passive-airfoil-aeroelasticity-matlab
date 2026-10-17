function report = run_cad_mechanisms(campaign)
%RUN_CAD_MECHANISMS Conditional original-solid A/B/C guide and load families.
% The slider is never detached from the imported monolithic solid. Guide
% coordinates constrain only its six mean rigid components at CAD faces 7/9.
root=fileparts(fileparts(mfilename('fullpath')));output=fullfile(campaign,'mechanisms');if ~isfolder(output),mkdir(output);end
rows=[];curves=[];
for level=1:3
 folder=fullfile(campaign,"tet10_level_"+level);s=load(fullfile(folder,'solution.mat'),'config','sensor','result');c=s.config;
 mesh=fem3d.elevate(fem3d.read_mesh(fullfile(root,string(c.discretization.meshes{level})),c.geometry.axisOrder,c.geometry.origin_m));
 system=fem3d.assemble(mesh,fem3d.material_reference(string(c.material.reference)));patch=cad.pressure_patch(mesh,c);mechanism=cad.slider_constraints(mesh,c,'B');
 fprintf('CAD guides level %d: six-coordinate compliance and gravity\n',level);
 gravity=fem3d.gravity(system,[0,0,-9.81]);F=[reshape((100*patch.unitForce)',[],1),reshape(gravity',[],1)];
 response=fem3d.constraint_response(system,F,s.result.bc,mechanism.rigidObservation);
 for scenario=1:3
  selected=[];if scenario==2,selected=1:6;elseif scenario==3,selected=2:6;end
  for loadScale=[0,.25,.5,1,.5,0]
   u=response.base(:,1)*loadScale;lambda=zeros(6,1);
   if ~isempty(selected),lambda(selected)=response.compliance(selected,selected)\(mechanism.rigidObservation(selected,:)*u);u=u-response.influence*lambda;end
   y=s.sensor.O*u;curves(end+1,:)=[level,scenario,loadScale*100,max(abs(y)),norm(lambda)]; %#ok<AGROW>
  end
  u=response.base(:,1);lambda=zeros(6,1);
  if ~isempty(selected),lambda(selected)=response.compliance(selected,selected)\(mechanism.rigidObservation(selected,:)*u);u=u-response.influence*lambda;end
  coordinates=mechanism.rigidObservation*u;y=s.sensor.O*u;
  rows(end+1,:)=[level,scenario,max(abs(y)),coordinates',max(response.energyResidual),norm(mechanism.rigidObservation(selected,:)*u,inf)]; %#ok<AGROW>
 end
 if level==3
  save(fullfile(output,'guide_compliance.mat'),'response','mechanism','c','-v7.3');
  gravityResult=struct('orientation',[0,0,-9.81],'maximumGravitySensor_m',max(abs(s.sensor.O*response.base(:,2))),'hypothesis','gravity along canonical height; reference orientation assumed');cad.write_json(fullfile(output,'gravity.json'),gravityResult);
 end
 clear system mesh response
end
report=array2table(rows,'VariableNames',{'level','scenario_ABC','maximum_TE_m','slider_x_m','slider_y_m','slider_z_m','slider_rx_rad','slider_ry_rad','slider_rz_rad','maximumEnergyResidual','constraintResidual'});
writetable(report,fullfile(output,'scenarios.csv'));writematrix(curves,fullfile(output,'pressure_loading_unloading.csv'));
cad.write_json(fullfile(output,'scope.json'),struct('assumed',true,'measured',false,'guide','candidate bar bounded by original CAD faces 7/9; mean chord translation permitted for C; all six blocked for B','loadCycles','elastic load/unload must coincide; friction separately verified, not measured hysteresis','connectivity','all original-solid connections retained','preload',0,'clearance','unknown; no invented finite guide clearance in A/B/C','experimentalValidation',false));disp(report);
end
