function op = action_operator(system)
%ACTION_OPERATOR Reference gradients for stable elementwise stiffness actions.
% Evaluating B' D B*u through strains avoids subtracting large K*u terms
% associated with local rigid translations. Same quadrature as assembled K.
T=system.mesh.tetra; X=system.mesh.nodes; ne=size(T,1); na=size(T,2);
if na==4
    reference={[-1,-1,-1;eye(3)]}; weights=1/6;
else
    a=(5+3*sqrt(5))/20;b=(5-sqrt(5))/20; points=[b,b,b;a,b,b;b,a,b;b,b,a];
    reference=cell(4,1); for k=1:4,[~,reference{k}]=fem3d.tet10_shape(points(k,:));end
    weights=ones(4,1)/24;
    edges=[1,2;2,3;3,1;1,4;2,4;3,4];
    for k=1:6
        mid=(X(T(:,edges(k,1)),:)+X(T(:,edges(k,2)),:))/2;
        assert(max(abs(X(T(:,k+4),:)-mid),[],'all')<1e-10,'fem3d:Action','Stable affine action requires straight Tet10 edges.');
    end
end
nq=numel(reference); G=cell(nq,1); for k=1:nq,G{k}=zeros(3,na,ne);end
w=zeros(ne,nq); rigid=zeros(3*na,6,ne);
for e=1:ne
    xe=X(T(e,:),:); xe=xe-xe(1,:);
    centered=xe-mean(xe,1);R=zeros(3*na,6);
    R(1:3:end,1)=1;R(2:3:end,2)=1;R(3:3:end,3)=1;
    for axis=1:3,omega=zeros(1,3);omega(axis)=1;motion=cross(repmat(omega,na,1),centered,2);R(:,axis+3)=reshape(motion',[],1);end
    R=R./vecnorm(R);[rigid(:,:,e),~]=qr(R,0);
    for k=1:nq
        J=xe'*reference{k}; grad=reference{k}/J;
        % Exact partition in the reference space; subtracting translations in
        % action() prevents absolute displacements amplifying roundoff.
        G{k}(:,:,e)=grad'; w(e,k)=abs(det(J))*weights(k);
    end
end
op=struct('gradient',{G},'rigid',rigid,'weights',w,'tetra',T,'materialId',system.materialId, ...
    'materials',system.materials,'nd',size(system.K,1),'na',na,'ne',ne);
end
