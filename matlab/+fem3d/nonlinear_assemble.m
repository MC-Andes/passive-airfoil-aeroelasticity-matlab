function result = nonlinear_assemble(system,u)
%NONLINEAR_ASSEMBLE Own finite-rotation solid internal force and tangent.
T=system.mesh.tetra;X=system.mesh.nodes;nd=3*size(X,1);ne=size(T,1);na=size(T,2);nl=3*na;
assert(numel(u)==nd && all(isfinite(u)),'fem3d:Displacement','Invalid nonlinear displacement.');
U=reshape(u,3,[])';force=zeros(nd,1);K=sparse(nd,nd);energy=0;minimumJ=Inf;
for start=1:400:ne
    ix=start:min(start+399,ne);I=zeros(nl^2,numel(ix));J=I;V=I;
    for j=1:numel(ix)
        e=ix(j);mat=system.materials(system.materialId(e));
        [f,k,w,detF]=fem3d.nonlinear_element(X(T(e,:),:),U(T(e,:),:),mat.D);
        dofs=reshape([3*T(e,:)-2;3*T(e,:)-1;3*T(e,:)],[],1);force(dofs)=force(dofs)+f;
        [ii,jj]=ndgrid(dofs,dofs);I(:,j)=ii(:);J(:,j)=jj(:);V(:,j)=k(:);energy=energy+w;minimumJ=min(minimumJ,detF);
    end
    K=K+sparse(I(:),J(:),V(:),nd,nd);
end
result=struct('force',force,'tangent',K,'strainEnergy',energy,'minimumDeformationJacobian',minimumJ, ...
    'constitutiveLaw','St Venant-Kirchhoff, total Lagrangian; elastic hypothesis');
end
