%RUN_FEM_CONVERGENCE Compare two supplied STEP volume meshes, same physics.
% This is a convergence diagnostic, not a claim of converged stresses.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'matlab'));
output=fullfile(root,'results','fem_reference'); if ~exist(output,'dir'), mkdir(output); end
[surface,~]=pae.surface_from_sections(fullfile(root,'data','reference'),4,14,"constant");
levels=["fem_reference","fem_refined"]; rows=zeros(2,8);
material=fem3d.material_reference("tough_pla_malim2023");
flow=struct('rho',1.225,'speed',10,'alpha',deg2rad(3));
for k=1:2
    mesh=fem3d.read_mesh(fullfile(root,'data',levels(k)),[1,3,2],surface.canonicalOriginCAD);
    system=fem3d.assemble(mesh,material);
    bc=struct('dofs',fem3d.dofs(find(mesh.nodes(:,1)<=.30*surface.referenceChord)));
    result=fem3d.one_way(system,surface,flow,bc);
    rows(k,:)=[size(mesh.nodes,1),size(mesh.tetra,1),system.totalMass, ...
        mesh.info.relative_volume_error,result.structure.maximumDisplacement*1e3, ...
        result.structure.stress.maximumVonMises/1e6,result.structure.strainEnergy, ...
        result.structure.relativeResidual];
end
convergence=array2table(rows,'VariableNames',{'nodes','tetra','mass_kg', ...
    'relative_volume_error','max_displacement_mm','max_vm_MPa','strain_energy_J','relative_residual'});
convergence.mesh=levels'; convergence=movevars(convergence,'mesh','Before',1);
writetable(convergence,fullfile(output,'convergence.csv'));
disp(convergence);
fprintf('Maximum displacement change relative to refined result: %.3g %%\n', ...
    100*abs(rows(2,5)-rows(1,5))/rows(2,5));
change=abs(rows(2,5)-rows(1,5))/rows(2,5);
if change>.05
    fprintf('NOT CONVERGED: displacement change exceeds the 5%% diagnostic tolerance.\n');
end
