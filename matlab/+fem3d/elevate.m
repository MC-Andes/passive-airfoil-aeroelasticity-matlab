function mesh = elevate(mesh)
%ELEVATE Conforming straight-edge Tet4 -> Tet10. No projection onto CAD.
assert(size(mesh.tetra,2)==4,'fem3d:Mesh','Elevation requires Tet4.');
T=mesh.tetra; ne=size(T,1); nn=size(mesh.nodes,1);
edges=[1,2;2,3;3,1;1,4;2,4;3,4]; allEdges=zeros(6*ne,2);
for k=1:6, allEdges((k-1)*ne+(1:ne),:)=sort(T(:,edges(k,:)),2); end
[E,~,g]=unique(allEdges,'rows');
mesh.nodes=[mesh.nodes;(mesh.nodes(E(:,1),:)+mesh.nodes(E(:,2),:))/2];
mesh.tetra=[T,reshape(nn+g,ne,6)];
if isfield(mesh,'cadFaces')
    F=mesh.cadFaces(:,1:3); nf=size(F,1); FE=[F(:,[1,2]);F(:,[2,3]);F(:,[3,1])];
    [ok,ids]=ismember(sort(FE,2),E,'rows'); assert(all(ok),'fem3d:Mesh','CAD face edge missing.');
    mesh.cadFaces=[F,reshape(nn+ids,nf,3)];
end
[mesh.faces,mesh.faceOwners]=fem3d.boundary(mesh.tetra,mesh.nodes);
mesh.elevation=struct('geometry','straight edges; no CAD reprojection','vertexCount',nn, ...
    'edges',E,'nodeOrder','vertices 1..4, edges 12 23 31 14 24 34');
end
