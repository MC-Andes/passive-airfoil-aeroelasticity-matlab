function system = assemble(mesh,materials)
%ASSEMBLE Sparse K/M from connectivity, coordinates and physical properties.
% mesh.nodes N x 3 [m]; mesh.tetra Ne x 4 (one-based).
% mesh.materialId optional Ne x 1, defaults to 1; materials is a struct array.
% A CAD volume label is not automatically a material ID.
X=mesh.nodes; T=mesh.tetra; nn=size(X,1); ne=size(T,1);
assert(size(X,2)==3 && all(isfinite(X(:))) && size(T,2)==4 && ne>0, ...
    'fem3d:Mesh','Require finite nodes and nonempty tetrahedra.');
assert(all(T(:)>=1 & T(:)<=nn & T(:)==fix(T(:))) && ...
    numel(unique(T(:)))==nn,'fem3d:Mesh','Invalid connectivity or orphan nodes.');
for m=1:numel(materials)
    D=materials(m).D;
    assert(isequal(size(D),[6,6]) && all(isfinite(D(:))) && ...
        norm(D-D','fro')<1e-10*norm(D,'fro'),'fem3d:Material','D must be finite and symmetric.');
    [~,p]=chol(D); assert(p==0,'fem3d:Material','D must be positive definite.');
    assert(isscalar(materials(m).rho) && isfinite(materials(m).rho) && ...
        materials(m).rho>0,'fem3d:Material','Require positive density.');
end
ids=ones(ne,1); if isfield(mesh,'materialId'), ids=mesh.materialId(:); end
assert(numel(ids)==ne && all(ids>=1 & ids<=numel(materials) & ids==fix(ids)), ...
    'fem3d:Material','Invalid material IDs.');
nd=3*nn; K=sparse(nd,nd); M=K; volumes=zeros(ne,1);
% Chunking avoids storing 144 triplets per element for an entire CAD mesh.
chunk=3000;
for start=1:chunk:ne
    last=min(start+chunk-1,ne); count=last-start+1;
    I=zeros(144,count); J=I; VK=I; VM=I;
    for e=start:last
        [ke,me,~,volumes(e)]=fem3d.tet4(X(T(e,:),:),materials(ids(e)).D,materials(ids(e)).rho);
        dofs=reshape([3*T(e,:)-2;3*T(e,:)-1;3*T(e,:)],1,12);
        [ii,jj]=ndgrid(dofs,dofs); c=e-start+1;
        I(:,c)=ii(:); J(:,c)=jj(:); VK(:,c)=ke(:); VM(:,c)=me(:);
    end
    K=K+sparse(I(:),J(:),VK(:),nd,nd);
    M=M+sparse(I(:),J(:),VM(:),nd,nd);
end
rho=arrayfun(@(m)m.rho,materials);
system=struct('K',(K+K')/2,'M',(M+M')/2,'mesh',mesh, ...
    'materials',materials,'materialId',ids,'volume',volumes, ...
    'totalMass',sum(volumes.*reshape(rho(ids),[],1)), ...
    'dofOrder','interleaved ux,uy,uz per structural node', ...
    'formulation','linear small-strain tetra4, consistent mass');
end
