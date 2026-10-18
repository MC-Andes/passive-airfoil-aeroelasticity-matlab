function [N,G] = tri_shape(r,count)
%TRI_SHAPE Tri3 or Tri6, vertices then edges [12 23 31].
L=[1-sum(r),r(:)']; dL=[-1,-1;1,0;0,1];
if count==3, N=L; G=dL; return; end
assert(count==6,'fem3d:Face','Require Tri3 or Tri6.');
N=zeros(1,6); G=zeros(6,2); N(1:3)=L.*(2*L-1); G(1:3,:)=(4*L'-1).*dL;
edges=[1,2;2,3;3,1];
for k=1:3
    a=edges(k,1); b=edges(k,2);
    N(k+3)=4*L(a)*L(b); G(k+3,:)=4*(L(a)*dL(b,:)+L(b)*dL(a,:));
end
end
