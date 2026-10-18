function bundle = build_reduced(solutionFolder,nx,ny,output)
%BUILD_REDUCED CAD material kinematics, verified modes and residual flexibility.
% K/M and a full one-way solve come from the own solid FEM. A residual vector
% is added only as a Galerkin basis vector, not as an arbitrary stiffness.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
s=load(fullfile(solutionFolder,'solution.mat'),'config','sensor','modes','result');config=s.config;
[~,name]=fileparts(solutionFolder);index=str2double(regexp(name,'\d+$','match','once'));
folder=fullfile(root,string(config.discretization.meshes{index}));
mesh=fem3d.read_mesh(folder,config.geometry.axisOrder,config.geometry.origin_m);
if size(mesh.tetra,2)==4,mesh=fem3d.elevate(mesh);end
fprintf('CAD reduction: material mean surface %dx%d\n',nx,ny);
surface=cad.surface(mesh,config,nx,ny);Phi=s.modes.Phi;
flow=struct('rho',config.flow.rho_kg_m3,'speed',config.flow.speed_m_s,'alpha',deg2rad(config.flow.alpha_deg), ...
    'referenceArea',surface.referenceArea,'referenceChord',surface.referenceChord, ...
    'referencePoint',[surface.referenceChord/4,0,mean(surface.nodes(1,:,3))]);
aero=fsi.vlm(surface.nodes,zeros(numel(surface.nodes)/3,3),flow);Fs=surface.H'*aero.nodalForces;
material=fem3d.material_reference(string(config.material.reference));system=fem3d.assemble(mesh,material);
fprintf('CAD reduction: full one-way static comparison\n');
full=fem3d.solve_static(system,Fs,s.result.bc);op=fem3d.action_operator(system);
Kr=diag(s.modes.omega.^2);modal=Phi*(Kr\(Phi'*reshape(Fs',[],1)));
yFull=s.sensor.O*full.u;yModal=s.sensor.O*modal;
errorBefore=norm(yModal-yFull,inf)/max(norm(yFull,inf),1e-15);
residual=full.u-Phi*(Phi'*system.M*full.u);mass=residual'*system.M*residual;
basis=Phi;if mass>1e-20*(full.u'*system.M*full.u),basis=[Phi,residual/sqrt(mass)];end
M=basis'*system.M*basis;action=fem3d.action(op,basis);K=basis'*action;K=(K+K')/2;
[W,D]=eig(K,M);[lambda,order]=sort(real(diag(D)));W=W(:,order);
for k=1:numel(lambda),W(:,k)=W(:,k)/sqrt(W(:,k)'*M*W(:,k));end
basis=basis*W;K=diag(lambda);M=eye(numel(lambda));O=s.sensor.O*basis;
B=surface.Oxyz*basis;zeta=.008;C=diag(2*zeta*sqrt(lambda));q=K\(basis'*reshape(Fs',[],1));
errorAfter=norm(O*q-yFull,inf)/max(norm(yFull,inf),1e-15);
assert(errorAfter<.02,'cad:Reduction','Residual-flexibility basis failed the full CAD one-way comparison.');
config.discretization.dampingRatio=zeta;config.provenance.runtimeOverrides='assumed zeta=0.008 from PLA literature';
model=struct('K',K,'M',M,'C',C,'PhiSurface',B,'X0',surface.nodes,'aero',@fsi.vlm, ...
    'observation',O,'frequencyHz',sqrt(lambda)/(2*pi),'zeta',zeta, ...
    'scope','CAD-derived Galerkin modes plus one residual-flexibility vector; assumed modal damping');
validation=struct('fullStaticSensorErrorBeforeCorrection',errorBefore,'fullStaticSensorErrorAfterCorrection',errorAfter, ...
    'fullEnergyResidual',full.energyResidual,'fullForceResidual',full.relativeResidual, ...
    'basisDimension',size(basis,2),'retainedEigenmodes',size(Phi,2), ...
    'staticLoad','rigid CAD VLM, same H/H-transpose material transfer','dampingSource','assumed 0.008 from PLA literature, not measured');
bundle=struct('model',model,'surface',surface,'config',config,'sensor',s.sensor,'validation',validation,'flow',flow);
if ~isfolder(output),mkdir(output);end
save(fullfile(output,'reduced_cad.mat'),'bundle','-v7.3');save(fullfile(output,'basis.mat'),'basis','-v7.3');
cad.write_json(fullfile(output,'reduction_validation.json'),validation);cad.provenance(root,config,output);
fprintf('CAD reduction: %d -> %d modes; full one-way sensor error %.3g -> %.3g\n',size(Phi,2),size(basis,2),errorBefore,errorAfter);
end
