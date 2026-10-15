function check_supports(mesh,bc)
%CHECK_SUPPORTS Every connected solid must have six rigid motions restrained.
% This check does not detect all internal mechanisms; Cholesky checks Kff too.
T=mesh.tetra; nn=size(mesh.nodes,1);
pairs=[1,2;1,3;1,4;2,3;2,4;3,4];
i=T(:,pairs(:,1)); j=T(:,pairs(:,2));
adj=sparse([i(:);j(:)],[j(:);i(:)],1,nn,nn);
component=conncomp(graph(adj,'upper'))';
for k=1:max(component)
    nodes=find(component==k); X=mesh.nodes(nodes,:);
    X=X-mean(X,1); X=X/max(max(vecnorm(X,2,2)),eps);
    R=zeros(3*numel(nodes),6);
    R(1:3:end,1)=1; R(2:3:end,2)=1; R(3:3:end,3)=1;
    for axis=1:3
        w=zeros(1,3); w(axis)=1;
        U=cross(repmat(w,numel(nodes),1),X,2);
        R(:,axis+3)=reshape(U',[],1);
    end
    globalDofs=fem3d.dofs(nodes); fixed=ismember(globalDofs,bc.dofs);
    assert(rank(R(fixed,:),1e-10)==6,'fem3d:Unconstrained', ...
        'Connected solid %d retains a rigid-body motion; specify physical supports.',k);
end
end
