function result = stress(system,u)
%STRESS Element centroid strains/stresses; no smoothing or allowable inference.
% von Mises is a stress summary, not a validated FDM failure criterion.
T=system.mesh.tetra; X=system.mesh.nodes; ne=size(T,1);
strain=zeros(ne,6); sigma=strain;
assert(numel(u)==3*size(X,1) && all(isfinite(u(:))),'fem3d:Displacement','Invalid displacement.');
[~,g10]=fem3d.tet10_shape([.25,.25,.25]);
for e=1:ne
    material=system.materials(system.materialId(e));
    if size(T,2)==4
        [~,~,B]=fem3d.tet4(X(T(e,:),:),material.D,material.rho);
    else
        xe=X(T(e,:),:); xe=xe-xe(1,:); B=fem3d.strain_matrix(g10/(xe'*g10));
    end
    dofs=reshape([3*T(e,:)-2;3*T(e,:)-1;3*T(e,:)],[],1);
    ue=reshape(u(dofs),3,[])';ue=ue-ue(1,:);
    strain(e,:)=(B*reshape(ue',[],1))'; sigma(e,:)=(material.D*strain(e,:)')';
end
vm=sqrt(.5*((sigma(:,1)-sigma(:,2)).^2+(sigma(:,2)-sigma(:,3)).^2+ ...
    (sigma(:,3)-sigma(:,1)).^2)+3*sum(sigma(:,4:6).^2,2));
result=struct('strain',strain,'sigma',sigma,'vonMises',vm, ...
    'maximumVonMises',max(vm),'location','material centroid; not a nodal maximum','componentOrder','xx yy zz xy yz xz (engineering strain shear)');
end
