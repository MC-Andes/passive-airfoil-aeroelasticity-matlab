function [Ke,Me,B,volume,info] = tet10(X,D,rho,order)
%TET10 Isoparametric quadratic solid; consistent mass, signed geometry checks.
% Straight edges: exact degree-2 stiffness and degree-4 mass. Curved edges:
% positive-weight Duffy quadrature, default 4^3; refine order to verify accuracy.
% B is evaluated at the material centroid, not a smoothed nodal stress.
arguments
    X (10,3) double {mustBeFinite}
    D (6,6) double {mustBeFinite}
    rho (1,1) double {mustBePositive,mustBeFinite}
    order (1,1) double {mustBeInteger,mustBePositive} = 4
end
X=X-X(1,:); % element-local origin avoids absolute CAD-coordinate cancellation
persistent Q W NN GG massReference
if isempty(Q)
    [x,w]=fem3d.gauss01(4); Q=zeros(64,3); W=zeros(64,1); NN=zeros(64,10); GG=zeros(10,3,64); k=0;
    for i=1:4, for j=1:4, for l=1:4
        k=k+1; Q(k,:)=[x(i),(1-x(i))*x(j),(1-x(i))*(1-x(j))*x(l)];
        W(k)=w(i)*w(j)*w(l)*(1-x(i))^2*(1-x(j));
        [NN(k,:),GG(:,:,k)]=fem3d.tet10_shape(Q(k,:));
    end, end, end
    massReference=NN'*(W.*NN);
end
edges=[1,2;2,3;3,1;1,4;2,4;3,4];
scale=max(vecnorm(X(2:4,:)-X(1,:),2,2));
straight=max(abs(X(5:10,:)-(X(edges(:,1),:)+X(edges(:,2),:))/2),[],'all')<1e-12*scale;
[~,gc]=fem3d.tet10_shape([.25,.25,.25]); J=X'*gc; d0=det(J);
assert(abs(d0)>1e-13*scale^3,'fem3d:Degenerate','Degenerate Tet10 centroid.');
B=fem3d.strain_matrix(gc/J); Ke=zeros(30); Me=Ke;
% Orientation may be negative (axis reflection); it must not change within an element.
controls=[0,0,0;1,0,0;0,1,0;0,0,1;.5,0,0;.5,.5,0;0,.5,0;0,0,.5;.5,0,.5;0,.5,.5];
minJ=abs(d0);
if ~straight
    for k=1:size(controls,1)
        [~,g]=fem3d.tet10_shape(controls(k,:)); d=det(X'*g);
        assert(d*d0>0 && abs(d)>1e-13*scale^3,'fem3d:Jacobian','Folded or degenerate Tet10 at a control point.');
        minJ=min(minJ,abs(d));
    end
end
if straight && nargin<4
    a=(5+3*sqrt(5))/20; b=(5-sqrt(5))/20;
    pts=[b,b,b;a,b,b;b,a,b;b,b,a];
    for k=1:4
        [~,g]=fem3d.tet10_shape(pts(k,:)); bk=fem3d.strain_matrix(g/J);
        Ke=Ke+abs(d0)/24*(bk'*D*bk);
    end
    volume=abs(d0)/6; Me=rho*abs(d0)*kron(massReference,eye(3));
    count=4;
else
    if order==4
        pts=Q; weights=W; shapes=NN; gradients=GG;
    else
        [x,w]=fem3d.gauss01(order); count=order^3;
        pts=zeros(count,3); weights=zeros(count,1); shapes=zeros(count,10); gradients=zeros(10,3,count); k=0;
        for i=1:order, for j=1:order, for l=1:order
            k=k+1; pts(k,:)=[x(i),(1-x(i))*x(j),(1-x(i))*(1-x(j))*x(l)];
            weights(k)=w(i)*w(j)*w(l)*(1-x(i))^2*(1-x(j));
            [shapes(k,:),gradients(:,:,k)]=fem3d.tet10_shape(pts(k,:));
        end, end, end
    end
    volume=0; count=size(pts,1);
    for k=1:count
        g=gradients(:,:,k); jk=X'*g; d=det(jk);
        assert(d*d0>0 && abs(d)>1e-13*scale^3,'fem3d:Jacobian','Folded or degenerate Tet10 at quadrature.');
        minJ=min(minJ,abs(d)); dv=weights(k)*abs(d); bk=fem3d.strain_matrix(g/jk);
        Ke=Ke+dv*(bk'*D*bk); Me=Me+rho*dv*kron(shapes(k,:)'*shapes(k,:),eye(3)); volume=volume+dv;
    end
end
Ke=(Ke+Ke')/2; Me=(Me+Me')/2;
info=struct('straightEdges',straight,'quadraturePoints',count,'minimumAbsJacobian',minJ, ...
    'orientation',sign(d0),'stressLocation','material centroid');
end
