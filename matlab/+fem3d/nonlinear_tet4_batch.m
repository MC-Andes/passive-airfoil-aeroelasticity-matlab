function [force,K,energy,minimumJ,fields] = nonlinear_tet4_batch(X,U,D)
%NONLINEAR_TET4_BATCH Same total-Lagrangian SVK as nonlinear_element, in pages.
% X,U: 4-by-3-by-Ne. Fields distinguish Green strain and Cauchy stress.
ne=size(X,3);assert(size(X,1)==4 && size(X,2)==3 && isequal(size(U),size(X)));
X=X-X(1,:,:);U=U-U(1,:,:);G=[-1,-1,-1;eye(3)];J=pagemtimes(permute(X,[2,1,3]),G);g=pagemrdivide(G,J);dv=abs(det3(J))/6;
F=eye(3)+pagemtimes(permute(U,[2,1,3]),g);jac=det3(F);minimumJ=min(jac);assert(all(jac>0),'fem3d:Inversion','Current configuration contains an inverted element.');
E=(pagemtimes(permute(F,[2,1,3]),F)-eye(3))/2;e=reshape(cat(1,E(1,1,:),E(2,2,:),E(3,3,:),2*E(1,2,:),2*E(2,3,:),2*E(1,3,:)),6,ne);s=D*e;
S=zeros(3,3,ne);S(1,1,:)=s(1,:);S(2,2,:)=s(2,:);S(3,3,:)=s(3,:);S(1,2,:)=s(4,:);S(2,1,:)=s(4,:);S(2,3,:)=s(5,:);S(3,2,:)=s(5,:);S(1,3,:)=s(6,:);S(3,1,:)=s(6,:);
B=zeros(6,12,ne);
for a=1:4,for r=1:3
 j=3*a-3+r;B(1,j,:)=g(a,1,:).*F(r,1,:);B(2,j,:)=g(a,2,:).*F(r,2,:);B(3,j,:)=g(a,3,:).*F(r,3,:);
 B(4,j,:)=g(a,2,:).*F(r,1,:)+g(a,1,:).*F(r,2,:);B(5,j,:)=g(a,3,:).*F(r,2,:)+g(a,2,:).*F(r,3,:);B(6,j,:)=g(a,3,:).*F(r,1,:)+g(a,1,:).*F(r,3,:);
end,end
force=reshape(pagemtimes(permute(B,[2,1,3]),reshape(s,6,1,ne)).*reshape(dv,1,1,ne),12,ne);
K=pagemtimes(permute(B,[2,1,3]),pagemtimes(D,B));geometric=pagemtimes(g,pagemtimes(S,permute(g,[2,1,3])));
for r=1:3,K(r:3:end,r:3:end,:)=K(r:3:end,r:3:end,:)+geometric;end
K=K.*reshape(dv,1,1,ne);energy=.5*sum(e.*s,1)'.*dv;
if nargout>4
 sigma=pagemtimes(F,pagemtimes(S,permute(F,[2,1,3])))./reshape(jac,1,1,ne);
 vm=sqrt(.5*((sigma(1,1,:)-sigma(2,2,:)).^2+(sigma(2,2,:)-sigma(3,3,:)).^2+(sigma(3,3,:)-sigma(1,1,:)).^2)+3*(sigma(1,2,:).^2+sigma(2,3,:).^2+sigma(1,3,:).^2));
 fields=struct('greenStrain',e','cauchyVonMises',reshape(vm,[],1),'deformationJacobian',jac);
end
end
function d=det3(A)
d=reshape(A(1,1,:).*(A(2,2,:).*A(3,3,:)-A(2,3,:).*A(3,2,:))-A(1,2,:).*(A(2,1,:).*A(3,3,:)-A(2,3,:).*A(3,1,:))+A(1,3,:).*(A(2,1,:).*A(3,2,:)-A(2,2,:).*A(3,1,:)),[],1);
end
