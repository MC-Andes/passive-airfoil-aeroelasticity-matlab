function [K,M,volume,minimumJacobian] = projected_operators(mesh,material,basis,order)
%PROJECTED_OPERATORS Isoparametric FEM Galerkin operators without full assembly.
% All elements and all positive Duffy quadrature points are integrated; this
% reduces storage, not the physics or quadrature. Batch arrays stay bounded.
[x,w]=fem3d.gauss01(order);nm=size(basis,2);K=zeros(nm);M=K;volume=0;minimumJacobian=inf;T=mesh.tetra;X=mesh.nodes;
for start=1:1000:size(T,1)
 ix=start:min(start+999,size(T,1));ne=numel(ix);nodes=T(ix,:);flat=nodes';flat=flat(:)';ids=reshape([3*flat-2;3*flat-1;3*flat],30,ne);
 xe=permute(reshape(X(nodes',:),10,ne,3),[1,3,2]);xe=xe-xe(1,:,:);
 ue=permute(reshape(basis(ids(:),:),3,10,ne,nm),[2,1,3,4]);ue=ue-ue(1,:,:,:);
 for i=1:order,for j=1:order,for k=1:order
  r=[x(i),(1-x(i))*x(j),(1-x(i))*(1-x(j))*x(k)];weight=w(i)*w(j)*w(k)*(1-x(i))^2*(1-x(j));[N,G]=fem3d.tet10_shape(r);
  J=pagemtimes(permute(xe,[2,1,3]),G);detJ=J(1,1,:).*(J(2,2,:).*J(3,3,:)-J(2,3,:).*J(3,2,:))-J(1,2,:).*(J(2,1,:).*J(3,3,:)-J(2,3,:).*J(3,1,:))+J(1,3,:).*(J(2,1,:).*J(3,2,:)-J(2,2,:).*J(3,1,:));
  assert(all(abs(detJ)>0,'all'),'fem3d:Jacobian','Zero projected jacobian.');minimumJacobian=min(minimumJacobian,min(abs(detJ),[],'all'));dv=reshape(abs(detJ)*weight,ne,1);
  grad=pagemrdivide(G,J);H=pagemtimes(permute(grad,[2,1,3]),ue);
  strain=cat(1,H(1,1,:,:),H(2,2,:,:),H(3,3,:,:),H(1,2,:,:)+H(2,1,:,:),H(2,3,:,:)+H(3,2,:,:),H(1,3,:,:)+H(3,1,:,:));
  E=reshape(permute(strain,[1,3,4,2]),6,ne*nm);stress=material.D*E;E=reshape(E,6*ne,nm);stress=reshape(stress,6*ne,nm);K=K+E'*(repelem(dv,6).*stress);
  % Absolute modal displacement is required for mass, unlike strain energy.
  absolute=reshape(basis(ids(:),:),3,10,ne,nm);vq=sum(absolute.*reshape(N,1,10,1,1),2);vq=reshape(vq,3*ne,nm);M=M+material.rho*vq'*(repelem(dv,3).*vq);volume=volume+sum(dv);
 end,end,end
end
K=(K+K')/2;M=(M+M')/2;
end
