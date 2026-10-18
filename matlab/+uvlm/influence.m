function A = influence(g,core)
%INFLUENCE Bound-ring no-penetration matrix.
np=size(g.ring,1); A=zeros(np);
if exist('aero_vortex_kernel','file')==3
 for first=1:64:np
  ix=first:min(first+63,np);v=zeros(np,3,numel(ix));
  for edge=1:4,v=v+uvlm.segment_influence(g.panels.control,g.ringVertices(g.ring(ix,edge),:),g.ringVertices(g.ring(ix,mod(edge,4)+1),:),core);end
  A(:,ix)=reshape(sum(v.*g.panels.normal,2),np,numel(ix));
 end
 assert(rcond(A)>1e-12,'uvlm:Lattice','Ill-conditioned bound-ring lattice.');return;
end
for j=1:np
    v=uvlm.ring_velocity(g.panels.control,g.ringVertices,g.ring(j,:),1,core);
    A(:,j)=sum(v.*g.panels.normal,2);
end
assert(rcond(A)>1e-12,'uvlm:Lattice','Ill-conditioned bound-ring lattice.');
end
