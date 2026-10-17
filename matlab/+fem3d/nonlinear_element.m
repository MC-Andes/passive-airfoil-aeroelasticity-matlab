function [force,tangent,energy,minimumJ] = nonlinear_element(X,U,D)
%NONLINEAR_ELEMENT Total-Lagrangian St Venant-Kirchhoff Tet4/Tet10.
% Finite rotations, Green strain and analytically consistent material/geometric
% tangent. Constitutive law is an explicit elastic hypothesis, not PLA damage.
na=size(X,1); assert(ismember(na,[4,10]) && isequal(size(U),size(X)),'fem3d:Element','Invalid solid state.');
X=X-X(1,:); U=U-U(1,:);
persistent points weights gradients
if isempty(points)
    [x,w]=fem3d.gauss01(4);points=zeros(64,3);weights=zeros(64,1);gradients=cell(64,1);k=0;
    for i=1:4,for j=1:4,for l=1:4
        k=k+1;points(k,:)=[x(i),(1-x(i))*x(j),(1-x(i))*(1-x(j))*x(l)];
        weights(k)=w(i)*w(j)*w(l)*(1-x(i))^2*(1-x(j));[~,gradients{k}]=fem3d.tet10_shape(points(k,:));
    end,end,end
end
if na==4, Gs={[-1,-1,-1;eye(3)]};W=1/6;else,Gs=gradients;W=weights;end
force=zeros(3*na,1);tangent=zeros(3*na);energy=0;minimumJ=Inf;
for k=1:numel(Gs)
    J=X'*Gs{k}; grad=Gs{k}/J; dV=abs(det(J))*W(k);
    F=eye(3)+U'*grad; j=det(F);minimumJ=min(minimumJ,j);
    assert(j>0,'fem3d:Inversion','Current configuration contains an inverted element.');
    E=(F'*F-eye(3))/2; e=[E(1,1);E(2,2);E(3,3);2*E(1,2);2*E(2,3);2*E(1,3)];
    stress=D*e;S=[stress(1),stress(4),stress(6);stress(4),stress(2),stress(5);stress(6),stress(5),stress(3)];
    B=zeros(6,3*na);
    for a=1:na
        g=grad(a,:);
        B(:,3*a-2:3*a)=[g(1)*F(:,1)';g(2)*F(:,2)';g(3)*F(:,3)'; ...
            g(2)*F(:,1)'+g(1)*F(:,2)';g(3)*F(:,2)'+g(2)*F(:,3)';g(3)*F(:,1)'+g(1)*F(:,3)'];
    end
    force=force+dV*(B'*stress);
    tangent=tangent+dV*(B'*D*B+kron(grad*S*grad',eye(3)));
    energy=energy+dV*.5*(e'*stress);
end
tangent=(tangent+tangent')/2;
end
