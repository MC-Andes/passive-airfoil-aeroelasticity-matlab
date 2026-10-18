function report = run_cad_gravity(solutionFolder)
%RUN_CAD_GRAVITY Orientation sensitivity, three independent full FEM loads.
root=fileparts(fileparts(mfilename('fullpath')));s=load(fullfile(solutionFolder,'solution.mat'),'config','sensor','result');c=s.config;
level=3;folder=fullfile(root,c.discretization.meshes{level});
out=fullfile(fileparts(solutionFolder),'gravity_orientation');if ~isfolder(out),mkdir(out);end;cad.provenance(root,c,out);
mesh=fem3d.elevate(fem3d.read_mesh(folder,c.geometry.axisOrder,c.geometry.origin_m));
system=fem3d.assemble(mesh,fem3d.material_reference(string(c.material.reference)));F=zeros(size(system.K,1),3);
for axis=1:3,g=zeros(1,3);g(axis)=-9.81;force=fem3d.gravity(system,g);F(:,axis)=reshape(force',[],1);end
r=fem3d.constraint_response(system,F,s.result.bc,sparse(0,size(system.K,1)));Y=s.sensor.O*r.base;
rows=zeros(3,5);for axis=1:3,U=reshape(r.base(:,axis),3,[])';rows(axis,:)=[axis,max(abs(Y(:,axis))),max(vecnorm(U,2,2)),r.energyResidual(axis),norm(sum(reshape(F(:,axis),3,[])',1))];end
report=array2table(rows,'VariableNames',{'gravity_axis_chord_span_height','maximum_TE_height_m','maximum_displacement_m','energy_residual','total_weight_N'});writetable(report,fullfile(out,'comparison.csv'));
scope=struct('gravityMagnitude_ms2',9.81,'primaryOrientationHypothesis','vertical span; gravity in canonical negative span direction','heightGravity','separate orientation sensitivity; cannot be silently combined with vertical-span primary case','zeroFlowReference','wind response is an increment about the zero-flow shape; finite prestress stiffening excluded by small-strain structural model','measured',false);cad.write_json(fullfile(out,'scope.json'),scope);save(fullfile(out,'solutions.mat'),'r','Y','report','scope','-v7.3');disp(report);
end
