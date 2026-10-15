function mesh = structured_mesh(grid)
%STRUCTURED_MESH Own mesher: split each logical hexahedron into six tets.
% grid(i,j,k,xyz), at least 2 nodes along each logical dimension.
sz=size(grid); assert(numel(sz)==4 && sz(4)==3 && all(sz(1:3)>=2), ...
    'fem3d:Mesh','grid must be ni x nj x nk x 3.');
mesh.nodes=reshape(grid,[],3); idx=reshape(1:prod(sz(1:3)),sz(1:3));
tetra=zeros(6*prod(sz(1:3)-1),4); row=0;
% Kuhn simplices all share the (000)-(111) body diagonal.
pattern=[1,2,4,8;1,4,3,8;1,3,7,8;1,7,5,8;1,5,6,8;1,6,2,8];
for k=1:sz(3)-1
    for j=1:sz(2)-1
        for i=1:sz(1)-1
            h=[idx(i,j,k),idx(i+1,j,k),idx(i,j+1,k),idx(i+1,j+1,k), ...
               idx(i,j,k+1),idx(i+1,j,k+1),idx(i,j+1,k+1),idx(i+1,j+1,k+1)];
            tetra(row+(1:6),:)=h(pattern); row=row+6;
        end
    end
end
mesh.tetra=tetra; mesh.materialId=ones(size(tetra,1),1);
[mesh.faces,mesh.faceOwners]=fem3d.boundary(tetra,mesh.nodes);
mesh.source="own structured tetrahedral mesher";
end
