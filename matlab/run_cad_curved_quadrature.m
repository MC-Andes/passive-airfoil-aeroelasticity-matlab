function report = run_cad_curved_quadrature(curvedFolder)
%RUN_CAD_CURVED_QUADRATURE Integrate actual curved-CAD solution fields at 4/6.
% This checks quadrature in the six computed curved modes and static field.
% It is not another full equilibrium solve or a bound over all FE vectors.
root=fileparts(fileparts(mfilename('fullpath')));
if nargin==0,curvedFolder=fullfile(root,'results/cad_original/meshes/curved_projection');end
a=load(fullfile(curvedFolder,'solution.mat'),'c','result','modes');
mesh=fem3d.read_mesh(curvedFolder,a.c.geometry.axisOrder,a.c.geometry.origin_m);
material=fem3d.material_reference(string(a.c.material.reference));
P=[a.modes.Phi,a.result.u];
% Normalize physical fields to keep their different units from hiding errors.
P=P./vecnorm(P);
operators=cell(2,1);orders=[4,6];
for i=1:2
    tic;[K,M,V,J]=fem3d.projected_operators(mesh,material,P,orders(i));
    operators{i}=struct('K',K,'M',M,'volume_m3',V,'minimumAbsoluteJacobian',J,'elapsed_s',toc);
    fprintf('Curved CAD quadrature order %d: %.1f s\n',orders(i),operators{i}.elapsed_s);
end
b=operators{2};dK=diag(1./sqrt(diag(b.K)));dM=diag(1./sqrt(diag(b.M)));
report=struct('orders',orders,'elements',size(mesh.tetra,1),'fields',size(P,2), ...
    'scaledStiffnessRelative',norm(dK*(operators{1}.K-b.K)*dK,'fro')/norm(dK*b.K*dK,'fro'), ...
    'scaledMassRelative',norm(dM*(operators{1}.M-b.M)*dM,'fro')/norm(dM*b.M*dM,'fro'), ...
    'fieldStrainEnergyRelative',abs(diag(operators{1}.K)./diag(b.K)-1), ...
    'fieldMassRelative',abs(diag(operators{1}.M)./diag(b.M)-1), ...
    'volumeRelative',abs(operators{1}.volume_m3/b.volume_m3-1), ...
    'minimumAbsoluteJacobian',[operators{1}.minimumAbsoluteJacobian,b.minimumAbsoluteJacobian], ...
    'scope','all actual curved CAD elements; six curved-solution modes plus curved static displacement, fixed fields at orders 4/6; quadrature diagnostic, not a full order-6 re-solve or certified operator bound');
report.targetRelative=.001;
report.passed=max([report.scaledStiffnessRelative;report.scaledMassRelative;report.fieldStrainEnergyRelative;report.fieldMassRelative;report.volumeRelative])<report.targetRelative;
cad.write_json(fullfile(curvedFolder,'quadrature_comparison.json'),report);
save(fullfile(curvedFolder,'quadrature_operators.mat'),'operators','report');disp(report);
end
