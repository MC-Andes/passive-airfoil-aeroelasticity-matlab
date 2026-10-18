function report = run_morphing_mechanics(casePath,order)
%RUN_MORPHING_MECHANICS Updated full CAD: shaft/servo tangent, modes, stress.
% A micrometre probe is not extrapolated into a finite design stroke.
root=fileparts(fileparts(mfilename('fullpath')));if nargin<1||isempty(casePath),casePath=fullfile(root,'cases/morphing_wing/case.json');end
c=morphing.load_case(root,casePath);if nargin<2,order=c.discretization.elementOrder;end
mesh=fem3d.read_mesh(fullfile(root,c.discretization.mesh),c.geometry.axisOrder,c.geometry.origin_m(:)');
assert(strcmp(mesh.info.source_sha256,c.geometry.sourceSha256),'morphing:Hash','Mesh differs from updated CAD.');
if order==10,mesh=fem3d.elevate(mesh);end
assert(abs(c.discretization.linearProbe_m)<=1e-4*c.geometry.chord_m,'morphing:LinearProbe','This tangent diagnostic cannot be used at design travel; use a verified finite-motion model.');
assert(ismember(order,[4,10]) && 3*size(mesh.nodes,1)<=c.discretization.maximumDofs,'morphing:Budget','Element order/resource contract failed.');
label=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS'));out=fullfile(root,'results/morphing_wing',['mechanics_',label]);mkdir(out);cad.provenance(root,c,out);
mat=struct('E',c.material.E_Pa,'nu',c.material.nu,'rho',c.material.rho_kg_m3,'D',fem3d.isotropic(c.material.E_Pa,c.material.nu));
fprintf('Updated CAD Tet%d assembly: %d DOFs\n',order,3*size(mesh.nodes,1));system=fem3d.assemble(mesh,mat);sensor=morphing.observers(mesh,c);
[bc,interface]=morphing.boundary_conditions(mesh,c,c.discretization.linearProbe_m);r=fem3d.solve_static(system,zeros(size(system.K,1),1),bc);
[bc0,~]=morphing.boundary_conditions(mesh,c,0);fprintf('Updated CAD locked-slider eigenmodes\n');modes=fem3d.modal(system,bc0,c.discretization.modalCount);
servo=sum(r.reaction(interface.driveNodes,:),1);shaft=sum(r.reaction(interface.mountNodes,:),1);boundaryWork=sum(reshape(r.displacement',[],1).*reshape(r.reaction',[],1));
assert(norm(servo+shaft)<1e-5*max(norm(servo),1) && abs(boundaryWork-2*r.strainEnergy)<1e-5*max(abs(boundaryWork),1e-12),'morphing:Balance','Driven boundary force/work balance failed.');
report=struct('output',out,'geometrySha256',c.geometry.sourceSha256,'order',order,'dofs',size(system.K,1),'elements',size(mesh.tetra,1), ...
 'mass_kg',system.totalMass,'command_m',c.discretization.linearProbe_m,'maxTE_m',max(abs(sensor.O*r.u)),'sensor_m',sensor.O*r.u, ...
 'driveForce_N',servo,'shaftForce_N',shaft,'energy_J',r.strainEnergy,'energyResidual',r.energyResidual,'assembledResidual',r.assembledRelativeResidual, ...
 'boundaryWorkDefect_J',boundaryWork-2*r.strainEnergy,'frequencyHz',modes.frequencyHz,'modalEnergyResidual',modes.eigenEnergyResidual, ...
 'maximumCentroidVonMises_Pa',max(r.stress.vonMises),'maximumCentroidStrain',max(abs(r.stress.strain),[],'all'), ...
 'scope','full updated CAD; linear tangent at zero-command, no aerodynamic forcing; PETG isotropic proxy and recommended insertion assumed; locked-slider free modes; stress diagnostic, no certified strength or mesh convergence');
cad.write_json(fullfile(out,'summary.json'),report);writetable(table(sensor.eta,sensor.points(:,1),sensor.O*r.u,'VariableNames',{'eta','x_m','w_m'}),fullfile(out,'sensors.csv'));writetable(table((1:numel(modes.frequencyHz))',modes.frequencyHz,modes.eigenEnergyResidual,'VariableNames',{'mode','frequency_Hz','energy_residual'}),fullfile(out,'modes.csv'));
save(fullfile(out,'solution.mat'),'c','r','modes','sensor','interface','report','-v7.3');
fig=fem3d.plot_result(system,r,'CAD actualizado PETG: prueba tangente del slider, no recorrido de diseño',1000);exportgraphics(fig,fullfile(out,'probe_deformation.png'),'Resolution',150);close(fig);disp(report);
end
