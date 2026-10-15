function [Ke,Me,B,volume] = tet4(X,D,rho)
%TET4 Own four-node solid FEM: exact constant-strain stiffness and mass.
% DOFs interleaved [ux1 uy1 uz1 ... ux4 uy4 uz4].
% Negative orientation is allowed; degenerate elements are rejected.
arguments
    X (4,3) double {mustBeFinite}
    D (6,6) double {mustBeFinite}
    rho (1,1) double {mustBePositive,mustBeFinite}
end
J=(X(2:4,:)-X(1,:))';
scale=max(vecnorm(J,2,1)); determinant=det(J);
assert(abs(determinant)>1e-13*scale^3,'fem3d:Degenerate','Degenerate tetrahedron.');
volume=abs(determinant)/6;
% Reference coordinates: N1=1-r-s-t, N2=r, N3=s, N4=t.
grad=[-1,-1,-1;1,0,0;0,1,0;0,0,1]/J;
B=zeros(6,12);
for a=1:4
    g=grad(a,:); idx=3*a-2:3*a;
    B(:,idx)=[g(1),0,0;0,g(2),0;0,0,g(3); ...
              g(2),g(1),0;0,g(3),g(2);g(3),0,g(1)];
end
Ke=volume*(B'*D*B); Ke=(Ke+Ke')/2;
% Integral Ni*Nj = V*(1+delta_ij)/20.
Me=rho*volume/20*kron(ones(4)+eye(4),eye(3));
end
