function result = stress(system,u)
%STRESS Constant element strains/stresses; no smoothing or allowable inference.
% von Mises is a stress summary, not a validated FDM failure criterion.
T=system.mesh.tetra; X=system.mesh.nodes; ne=size(T,1);
strain=zeros(ne,6); sigma=strain;
assert(numel(u)==3*size(X,1) && all(isfinite(u(:))),'fem3d:Displacement','Invalid displacement.');
for e=1:ne
    material=system.materials(system.materialId(e));
    [~,~,B]=fem3d.tet4(X(T(e,:),:),material.D,material.rho);
    dofs=reshape([3*T(e,:)-2;3*T(e,:)-1;3*T(e,:)],[],1);
    strain(e,:)=(B*u(dofs))'; sigma(e,:)=(material.D*strain(e,:)')';
end
vm=sqrt(.5*((sigma(:,1)-sigma(:,2)).^2+(sigma(:,2)-sigma(:,3)).^2+ ...
    (sigma(:,3)-sigma(:,1)).^2)+3*sum(sigma(:,4:6).^2,2));
result=struct('strain',strain,'sigma',sigma,'vonMises',vm, ...
    'maximumVonMises',max(vm),'componentOrder','xx yy zz xy yz xz (engineering strain shear)');
end
