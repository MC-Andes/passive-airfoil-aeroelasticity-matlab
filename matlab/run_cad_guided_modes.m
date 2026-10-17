function report = run_cad_guided_modes(campaign,level,count)
%RUN_CAD_GUIDED_MODES Full constrained CAD modes, no plate replacement.
if nargin<2,level=3;end;if nargin<3,count=3;end
root=fileparts(fileparts(mfilename('fullpath')));s=load(fullfile(campaign,"tet10_level_"+level,'solution.mat'),'config','result');c=s.config;mesh=fem3d.elevate(fem3d.read_mesh(fullfile(root,string(c.discretization.meshes{level})),c.geometry.axisOrder,c.geometry.origin_m));system=fem3d.assemble(mesh,fem3d.material_reference(string(c.material.reference)));rows=[];
for scenario=['B','C']
 mechanism=cad.slider_constraints(mesh,c,scenario);fprintf('CAD full guided modes scenario %s\n',scenario);modes=fem3d.modal_constrained(system,s.result.bc,mechanism.C,count);save(fullfile(campaign,'mechanisms',"modes_"+scenario+"_level_"+level+'.mat'),'modes','mechanism','-v7.3');
 for k=1:count,rows(end+1,:)=[double(scenario)-64,k,modes.frequencyHz(k),modes.eigenEnergyResidual(k),modes.constraintResidual(k)];end %#ok<AGROW>
end
report=array2table(rows,'VariableNames',{'scenario_ABC','mode','frequency_Hz','energy_residual','constraint_residual'});writetable(report,fullfile(campaign,'mechanisms',"guided_modes_level_"+level+'.csv'));disp(report);
end
