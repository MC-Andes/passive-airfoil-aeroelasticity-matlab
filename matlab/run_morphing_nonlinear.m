function report = run_morphing_nonlinear(command,casePath)
%RUN_MORPHING_NONLINEAR Actual updated CAD with prescribed finite-motion BC.
% A conditional SVK calculation; it does not establish PETG failure strength.
root=fileparts(fileparts(mfilename('fullpath')));if nargin<2||isempty(casePath),casePath=fullfile(root,'cases/morphing_wing/case.json');end
c=morphing.load_case(root,casePath);if nargin<1,command=c.discretization.smallNonlinearProbe_m;end
mesh=fem3d.read_mesh(fullfile(root,c.discretization.mesh),c.geometry.axisOrder,c.geometry.origin_m(:)');mesh.nonlinearBatch=true;
mat=struct('E',c.material.E_Pa,'nu',c.material.nu,'rho',c.material.rho_kg_m3,'D',fem3d.isotropic(c.material.E_Pa,c.material.nu));
label=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS'));out=fullfile(root,'results/morphing_wing',['nonlinear_',label]);mkdir(out);cad.provenance(root,c,out);
fprintf('Updated CAD nonlinear: prescribed slider %g mm\n',1000*command);system=fem3d.assemble(mesh,mat);[bc,interface]=morphing.boundary_conditions(mesh,c,command);sensor=morphing.observers(mesh,c);
r=fem3d.solve_nonlinear(system,zeros(size(system.K,1),1),bc,struct('loadSteps',4,'maximumIterations',25,'relativeTolerance',1e-7));
T=mesh.tetra;maximumGreen=0;vm=zeros(size(T,1),1);J=vm;
for first=1:2000:size(T,1)
 ids=first:min(first+1999,size(T,1));nodes=T(ids,:);X=permute(reshape(mesh.nodes(nodes',:),4,numel(ids),3),[1,3,2]);U=permute(reshape(r.displacement(nodes',:),4,numel(ids),3),[1,3,2]);
 [~,~,~,~,f]=fem3d.nonlinear_tet4_batch(X,U,mat.D);maximumGreen=max(maximumGreen,max(abs(f.greenStrain),[],'all'));vm(ids)=f.cauchyVonMises;J(ids)=f.deformationJacobian;
end
force=reshape(r.state.force,3,[])';servo=sum(force(interface.driveNodes,:),1);shaft=sum(force(interface.mountNodes,:),1);
report=struct('output',out,'command_m',command,'maxTE_m',max(abs(sensor.O*r.u)),'maxGreenStrain',maximumGreen,'minimumDeformationJacobian',min(J), ...
 'maximumCauchyVonMises_Pa',max(vm),'driveForce_N',servo,'shaftForce_N',shaft,'strainEnergy_J',r.state.strainEnergy,'iterations',size(r.history,1),'converged',r.converged, ...
 'geometrySha256',c.geometry.sourceSha256,'scope','actual updated Tet4 CAD, conditional SVK and rigid driven face; no flow, no mesh convergence, no calibrated PETG plasticity/damage; neither a strength qualification nor permission to linearly extrapolate design travel');
assert(norm(servo+shaft)<1e-4*max(norm(servo),1),'morphing:Balance','Nonlinear support/servo balance failed.');
cad.write_json(fullfile(out,'summary.json'),report);writematrix(r.history,fullfile(out,'newton_history.csv'));save(fullfile(out,'solution.mat'),'r','sensor','interface','report','c','vm','J','-v7.3');disp(report);
end
