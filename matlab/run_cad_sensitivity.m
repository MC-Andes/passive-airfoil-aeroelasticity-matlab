function report = run_cad_sensitivity(solutionFolder)
%RUN_CAD_SENSITIVITY Isolated support and height/thickness geometric variants.
% Original CAD is unchanged. Geometry variants are explicit affine maps of
% its cap-derived mesh and are never substituted for the principal case.
root=fileparts(fileparts(mfilename('fullpath')));s=load(fullfile(solutionFolder,'solution.mat'),'config');c=s.config;[~,name]=fileparts(solutionFolder);level=str2double(regexp(name,'\d+$','match','once'));
original=fem3d.elevate(fem3d.read_mesh(fullfile(root,string(c.discretization.meshes{level})),c.geometry.axisOrder,c.geometry.origin_m));rows=[];variants=[1,.30;1,.25;1,.35;.95,.30;1.05,.30];
output=fullfile(fileparts(solutionFolder),'sensitivity');if ~isfolder(output),mkdir(output);end
for k=1:size(variants,1)
 mesh=original;z0=mean([min(mesh.nodes(:,3)),max(mesh.nodes(:,3))]);mesh.nodes(:,3)=z0+variants(k,1)*(mesh.nodes(:,3)-z0);config=c;config.supports.x_over_c=variants(k,2);
 system=fem3d.assemble(mesh,fem3d.material_reference(string(c.material.reference)));p=cad.pressure_patch(mesh,config);sensor=cad.sensors(mesh,config);bc=struct('dofs',fem3d.dofs(find(mesh.nodes(:,1)<=variants(k,2)*c.geometry.chord_m)));r=fem3d.solve_static(system,100*p.unitForce,bc);m=fem3d.modal(system,bc,1);
 rows(end+1,:)=[k,variants(k,:),max(abs(sensor.O*r.u)),system.totalMass,m.frequencyHz,r.energyResidual,m.eigenEnergyResidual]; %#ok<AGROW>
 fprintf('CAD sensitivity %d: height scale %.3g, clamp %.3g, TE %.5g mm, f %.4g Hz\n',k,variants(k,1),variants(k,2),1e3*rows(end,4),m.frequencyHz);
end
report=array2table(rows,'VariableNames',{'variant','canonical_height_scale','clamp_x_over_c','maximum_TE_m','mass_kg','frequency_Hz','static_energy_residual','modal_energy_residual'});writetable(report,fullfile(output,'support_thickness.csv'));
cad.write_json(fullfile(output,'scope.json'),struct('originalCADAltered',false,'referenceMesh',string(c.discretization.meshes{level}),'geometryVariation','explicit uniform canonical-height scale about fixed midheight; skin/ribs/cavity heights vary together; not a claim of isolated wall-gauge tolerance','thicknessRange','engineering ±5%; not a measured manufacturing distribution','supports','clamp boundary at 25/30/35% chord, no support springs inferred','meshError','compare baseline with the independent three-level campaign; these are sensitivity variants, not new validated geometries','experimentalValidation',false));
end
