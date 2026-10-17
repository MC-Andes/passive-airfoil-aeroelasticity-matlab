function [K,M,volume,minimumRatio] = tet10_batch(X,D,rho,order)
%TET10_BATCH Bounded-memory vectorized isoparametric Tet10 integration.
% X is 10 x 3 x Ne; each page has the same material. Same positive Duffy rule
% and extra Jacobian controls as tet10; intended for curved CAD campaigns.
if nargin<4,order=4;end
ne=size(X,3);X=X-X(1,:,:);[x,w]=fem3d.gauss01(order);K=zeros(30,30,ne);M=K;volume=zeros(ne,1);minimumRatio=inf;
[~,G]=fem3d.tet10_shape([.25,.25,.25]);centroid=determinant(pagemtimes(permute(X,[2,1,3]),G));assert(all(centroid~=0));
controls=[0,0,0;1,0,0;0,1,0;0,0,1;.5,0,0;.5,.5,0;0,.5,0;0,0,.5;.5,0,.5;0,.5,.5];
for k=1:size(controls,1),[~,G]=fem3d.tet10_shape(controls(k,:));check(pagemtimes(permute(X,[2,1,3]),G));end
for i=1:order,for j=1:order,for k=1:order
 r=[x(i),(1-x(i))*x(j),(1-x(i))*(1-x(j))*x(k)];weight=w(i)*w(j)*w(k)*(1-x(i))^2*(1-x(j));[N,G]=fem3d.tet10_shape(r);J=pagemtimes(permute(X,[2,1,3]),G);d=check(J);dv=abs(d)*weight;
 grad=pagemrdivide(G,J);B=zeros(6,30,ne);
 B(1,1:3:end,:)=permute(grad(:,1,:),[2,1,3]);B(2,2:3:end,:)=permute(grad(:,2,:),[2,1,3]);B(3,3:3:end,:)=permute(grad(:,3,:),[2,1,3]);
 B(4,1:3:end,:)=permute(grad(:,2,:),[2,1,3]);B(4,2:3:end,:)=permute(grad(:,1,:),[2,1,3]);B(5,2:3:end,:)=permute(grad(:,3,:),[2,1,3]);B(5,3:3:end,:)=permute(grad(:,2,:),[2,1,3]);B(6,1:3:end,:)=permute(grad(:,3,:),[2,1,3]);B(6,3:3:end,:)=permute(grad(:,1,:),[2,1,3]);
 K=K+pagemtimes(permute(B,[2,1,3]),pagemtimes(D,B)).*reshape(dv,1,1,ne);M=M+rho*kron(N'*N,eye(3)).*reshape(dv,1,1,ne);volume=volume+dv;
end,end,end
 function d=check(J)
 d=determinant(J);assert(all(d.*centroid>0),'fem3d:Jacobian','Folded curved Tet10.');minimumRatio=min(minimumRatio,min(abs(d./centroid)));
 end
end
function d=determinant(J)
d=reshape(J(1,1,:).*(J(2,2,:).*J(3,3,:)-J(2,3,:).*J(3,2,:))-J(1,2,:).*(J(2,1,:).*J(3,3,:)-J(2,3,:).*J(3,1,:))+J(1,3,:).*(J(2,1,:).*J(3,2,:)-J(2,2,:).*J(3,1,:)),[],1);
end
