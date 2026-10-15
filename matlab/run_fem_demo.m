%RUN_FEM_DEMO From-scratch solid FEM + existing own steady VLM; one-way.
% All numerical solvers use MATLAB base. Parametric solid NACA surrogate.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'matlab'));
output=fullfile(root,'results','fem_demo'); if ~exist(output,'dir'), mkdir(output); end
material=fem3d.material_reference("tough_pla_malim2023");
mesh=fem3d.airfoil_mesh(.247,.36,.12,[30,16,4]);
system=fem3d.assemble(mesh,material);
% Assumed tunnel fixture: the anterior 30% chord is clamped throughout span.
fixed=find(mesh.nodes(:,1)<=.3*mesh.chord);
bc=struct('dofs',fem3d.dofs(fixed),'description',"ASSUMED clamped anterior 30% chord");
p=struct('spanStations',[-.18;.18],'leadingEdge',[0,-.18,0;0,.18,0], ...
    'chord',[.247;.247],'camberFcn',@(xi,eta) zeros(size(xi)));
surface=pae.make_surface(p,6,18);
flow=struct('rho',1.225,'speed',10,'alpha',deg2rad(3));
result=fem3d.one_way(system,surface,flow,bc);
modes=fem3d.modal(system,bc,6);
% A small structural free vibration in the first assembled mode, no fluid feedback.
u0=modes.Phi(:,1); u0=u0*(1e-5/max(abs(u0)));
time=linspace(0,3/modes.frequencyHz(1),401);
transient=fem3d.newmark(system,bc,time,@(t) zeros(size(system.K,1),1),u0);
summary=table(string(material.name),size(mesh.nodes,1),size(mesh.tetra,1),system.totalMass, ...
    flow.speed,rad2deg(flow.alpha),result.aero.CL,result.aero.lift, ...
    result.structure.maximumDisplacement*1e3,result.structure.stress.maximumVonMises/1e6, ...
    result.structure.relativeResidual,modes.frequencyHz(1), ...
    'VariableNames',{'material','nodes','tetra','mass_kg','U_ms','alpha_deg','CL', ...
    'lift_N','max_displacement_mm','max_vm_MPa','relative_residual','f1_Hz'});
writetable(summary,fullfile(output,'summary.csv'));
writetable(table((1:6)',modes.frequencyHz,modes.eigenResidual, ...
    'VariableNames',{'mode','frequency_Hz','residual'}),fullfile(output,'modes.csv'));
writetable(table(time',transient.energy','VariableNames',{'time_s','energy_J'}),fullfile(output,'energy.csv'));
fig=fem3d.plot_result(system,result.structure, ...
    "FEM propio | NACA macizo | Tough PLA | apoyo supuesto | cargas desacopladas",20);
exportgraphics(fig,fullfile(output,'fem_demo.png'),'Resolution',160); close(fig);
save(fullfile(output,'study.mat'),'system','surface','material','result','modes','transient','summary','-v7.3');
disp(summary);
