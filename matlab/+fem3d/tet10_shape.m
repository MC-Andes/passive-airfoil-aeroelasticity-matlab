function [N,G] = tet10_shape(r)
%TET10_SHAPE Quadratic tetrahedron, local r=[r,s,t].
% Nodes: vertices 1..4; edges [12 23 31 14 24 34]. G=dN/d[r,s,t].
L=[1-sum(r),r(:)']; dL=[-1,-1,-1;eye(3)];
N=zeros(1,10); G=zeros(10,3); N(1:4)=L.*(2*L-1);
G(1:4,:)=(4*L'-1).*dL;
edges=[1,2;2,3;3,1;1,4;2,4;3,4];
for k=1:6
    a=edges(k,1); b=edges(k,2);
    N(k+4)=4*L(a)*L(b); G(k+4,:)=4*(L(a)*dL(b,:)+L(b)*dL(a,:));
end
end
