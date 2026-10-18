function q = surface_quadrature(mesh,faces,order)
%SURFACE_QUADRATURE Consistent kinematics on explicitly supplied physical faces.
% q.H maps nodes to integration points. q.H' is its work-conjugate force map.
% q.weights are actual curved area weights [m^2]; normals are outward if faces are.
arguments
    mesh struct
    faces double
    order (1,1) double {mustBeInteger,mustBePositive} = 4
end
nf=size(faces,1); nc=size(faces,2); nn=size(mesh.nodes,1);
assert(nf>0 && ismember(nc,[3,6]) && all(faces(:)>=1 & faces(:)<=nn & faces(:)==fix(faces(:))), ...
    'fem3d:Face','Invalid/nonempty Tri3 or Tri6 faces required.');
[x,w]=fem3d.gauss01(order); nq=order^2; shape=zeros(nq,nc); grad=zeros(nc,2,nq); wt=zeros(nq,1); k=0;
for i=1:order, for j=1:order
    k=k+1; [shape(k,:),grad(:,:,k)]=fem3d.tri_shape([x(i),(1-x(i))*x(j)],nc);
    wt(k)=w(i)*w(j)*(1-x(i));
end, end
I=zeros(nf*nq,nc); J=I; V=I; normal=zeros(nf*nq,3); weights=zeros(nf*nq,1);
for f=1:nf
    rows=(f-1)*nq+(1:nq); X=mesh.nodes(faces(f,:),:);
    I(rows,:)=repmat(rows',1,nc); J(rows,:)=repmat(faces(f,:),nq,1); V(rows,:)=shape;
    for k=1:nq
        tangent=X'*grad(:,:,k); n=cross(tangent(:,1),tangent(:,2))'; jac=norm(n);
        assert(jac>0,'fem3d:Face','Degenerate face at integration point.');
        normal(rows(k),:)=n/jac; weights(rows(k))=wt(k)*jac;
    end
end
H=sparse(I(:),J(:),V(:),nf*nq,nn);
q=struct('H',H,'points',H*mesh.nodes,'normals',normal,'weights',weights, ...
    'faceIndex',repelem((1:nf)',nq),'area',sum(weights),'quadratureOrder',order);
end
