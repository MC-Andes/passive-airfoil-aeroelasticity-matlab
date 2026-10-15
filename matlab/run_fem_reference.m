%RUN_FEM_REFERENCE Own solid FEM on supplied STEP volume; one-way VLM loads.
% Hypothetical homogeneous bonded solid with assumed supports. NOT the actual
% slider mechanism: STEP is one BREP and does not supply joints or constraints.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'matlab'));
output=fullfile(root,'results','fem_reference'); if ~exist(output,'dir'), mkdir(output); end
[surface,~]=pae.surface_from_sections(fullfile(root,'data','reference'),4,14,"constant");
mesh=fem3d.read_mesh(fullfile(root,'data','fem_reference'),[1,3,2],surface.canonicalOriginCAD);
material=fem3d.material_reference("tough_pla_malim2023");
system=fem3d.assemble(mesh,material);
% Assumed support, explicit and editable. It is not identified from the CAD.
fixed=find(mesh.nodes(:,1)<=.30*surface.referenceChord);
bc=struct('dofs',fem3d.dofs(fixed),'description',"ASSUMED clamped anterior 30% chord; one bonded solid");
flow=struct('rho',1.225,'speed',10,'alpha',deg2rad(3));
result=fem3d.one_way(system,surface,flow,bc);
result.meshConvergenceValidated=false;
result.interpretation="PRELIMINARY: assumed material/supports; STEP mesh displacement not converged";
summary=table(string(material.name),size(mesh.nodes,1),size(mesh.tetra,1),system.totalMass, ...
    flow.speed,rad2deg(flow.alpha),result.aero.CL,result.aero.lift, ...
    result.structure.maximumDisplacement*1e3,result.structure.stress.maximumVonMises/1e6, ...
    result.structure.relativeResidual,norm(result.structure.forceBalance), ...
    norm(result.structure.momentBalance),result.transferInfo.maximumPatchRadius, ...
    'VariableNames',{'material','nodes','tetra','mass_kg','U_ms','alpha_deg','CL', ...
    'lift_N','max_displacement_mm','max_vm_MPa','relative_residual','force_balance_N', ...
    'moment_balance_Nm','max_transfer_radius_m'});
writetable(summary,fullfile(output,'summary.csv'));
writematrix([mesh.nodes,result.structure.displacement],fullfile(output,'displacement_xyz_m.csv'));
fig=fem3d.plot_result(system,result.structure, ...
    "STEP: resultado sin convergencia de malla | material y apoyo supuestos | sin cinemática del slider",20);
exportgraphics(fig,fullfile(output,'fem_reference.png'),'Resolution',160); close(fig);
save(fullfile(output,'study.mat'),'system','surface','material','result','summary','-v7.3');
disp(summary);
