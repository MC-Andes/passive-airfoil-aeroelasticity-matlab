function report = run_cad_full_reduced(solutionFolder)
%RUN_CAD_FULL_REDUCED Actual CAD full vs reduced free-mode response and energy.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'matlab'));
s=load(fullfile(solutionFolder,'solution.mat'),'config','sensor','modes','result');[~,name]=fileparts(solutionFolder);ix=str2double(regexp(name,'\d+$','match','once'));
mesh=fem3d.read_mesh(fullfile(root,string(s.config.discretization.meshes{ix})),s.config.geometry.axisOrder,s.config.geometry.origin_m);if size(mesh.tetra,2)==4,mesh=fem3d.elevate(mesh);end
system=fem3d.assemble(mesh,fem3d.material_reference(string(s.config.material.reference)));omega=s.modes.omega(1);zeta=.008;
u0=s.modes.Phi(:,1);u0=u0*1e-5/max(abs(s.sensor.O*u0));time=linspace(0,4/s.modes.frequencyHz(1),401);alpha=2*zeta*omega;
full=fem3d.newmark_observed(system,s.result.bc,time,@(t)zeros(size(u0)),s.sensor.O,u0,zeros(size(u0)),[alpha,0]);
model=struct('K',omega^2,'M',1,'C',alpha,'B',0,'O',s.sensor.O*u0);reduced=ael.transient(model,time,@(t)0,1,0);
error=norm(full.Y-reduced.Y,inf)/max(norm(reduced.Y,inf),1e-15);
report=struct('maximumRowSumRelativeError',error,'maximumPointwiseRelativeError',max(abs(full.Y-reduced.Y),[],'all')/1e-5, ...
    'fullRelativeEnergyDefect',full.relativeEnergyBalance,'fullForceResidual_N',full.maximumForceResidual, ...
    'steps',400,'periods',4,'dampingRatioAssumed',zeta,'scope','actual full CAD and single-mode exact-subspace structural comparison; no flow forcing');
assert(report.maximumPointwiseRelativeError<.02 && full.relativeEnergyBalance<.02,'cad:Reduction','Full/reduced CAD dynamic comparison failed.');
out=fullfile(solutionFolder,'full_reduced');if ~isfolder(out),mkdir(out);end
cad.write_json(fullfile(out,'comparison.json'),report);save(fullfile(out,'response.mat'),'full','reduced','report','-v7.3');
writematrix([time',full.Y',reduced.Y'],fullfile(out,'sensors.csv'));disp(report);
end
