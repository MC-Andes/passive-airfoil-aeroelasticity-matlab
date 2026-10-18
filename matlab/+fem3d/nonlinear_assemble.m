function result = nonlinear_assemble(system,u)
%NONLINEAR_ASSEMBLE Own finite-rotation solid internal force and tangent.
T=system.mesh.tetra;X=system.mesh.nodes;nd=3*size(X,1);ne=size(T,1);na=size(T,2);nl=3*na;
assert(numel(u)==nd && all(isfinite(u)),'fem3d:Displacement','Invalid nonlinear displacement.');
U=reshape(u,3,[])';force=zeros(nd,1);K=sparse(nd,nd);energy=0;minimumJ=Inf;
batch=na==4 && isfield(system.mesh,'nonlinearBatch') && system.mesh.nonlinearBatch;
chunk=400;if batch,chunk=2000;end
for start=1:chunk:ne
    ix=start:min(start+chunk-1,ne);I=zeros(nl^2,numel(ix));J=I;V=I;
    if batch
        nodes=T(ix,:);flat=nodes';flat=flat(:)';ids=reshape([3*flat-2;3*flat-1;3*flat],12,[]);
        xe=permute(reshape(X(nodes',:),4,numel(ix),3),[1,3,2]);ue=permute(reshape(U(nodes',:),4,numel(ix),3),[1,3,2]);
        for material=unique(system.materialId(ix))'
            use=find(system.materialId(ix)==material);[f,k,w,detF]=fem3d.nonlinear_tet4_batch(xe(:,:,use),ue(:,:,use),system.materials(material).D);
            force=force+accumarray(reshape(ids(:,use),[],1),f(:),[nd,1]);V(:,use)=reshape(k,nl^2,[]);energy=energy+sum(w);minimumJ=min(minimumJ,detF);
        end
        for j=1:numel(ix),[ii,jj]=ndgrid(ids(:,j),ids(:,j));I(:,j)=ii(:);J(:,j)=jj(:);end
    else
    for j=1:numel(ix)
        e=ix(j);mat=system.materials(system.materialId(e));
        [f,k,w,detF]=fem3d.nonlinear_element(X(T(e,:),:),U(T(e,:),:),mat.D);
        dofs=reshape([3*T(e,:)-2;3*T(e,:)-1;3*T(e,:)],[],1);force(dofs)=force(dofs)+f;
        [ii,jj]=ndgrid(dofs,dofs);I(:,j)=ii(:);J(:,j)=jj(:);V(:,j)=k(:);energy=energy+w;minimumJ=min(minimumJ,detF);
    end
    end
    K=K+sparse(I(:),J(:),V(:),nd,nd);
end
result=struct('force',force,'tangent',K,'strainEnergy',energy,'minimumDeformationJacobian',minimumJ, ...
    'constitutiveLaw','St Venant-Kirchhoff, total Lagrangian; elastic hypothesis');
end
