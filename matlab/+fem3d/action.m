function force = action(op,U)
%ACTION Stable multi-vector internal-force evaluation of linear elasticity.
% Identical linear operator in exact arithmetic; local translation removed
% before differentiation. No stiffness regularization or artificial springs.
assert(size(U,1)==op.nd,'fem3d:Displacement','Invalid action vector.');
nm=size(U,2); force=zeros(op.nd,nm); na=op.na;
for start=1:2500:op.ne
    ix=start:min(start+2499,op.ne); ne=numel(ix); T=op.tetra(ix,:);
    flat=T'; flat=flat(:)'; ids=reshape([3*flat-2;3*flat-1;3*flat],3*na,ne);
    rigid=op.rigid(:,:,ix); ue=permute(reshape(U(ids(:),:),3*na,ne,nm),[1,3,2]);
    ue=ue-pagemtimes(rigid,pagemtimes(permute(rigid,[2,1,3]),ue));
    ue=permute(reshape(permute(ue,[1,3,2]),3,na,ne,nm),[2,1,3,4]);
    local=zeros(na,3,ne,nm); D=zeros(6,6,ne);
    for material=1:numel(op.materials),D(:,:,op.materialId(ix)==material)=repmat(op.materials(material).D,1,1,sum(op.materialId(ix)==material));end
    for q=1:numel(op.gradient)
        G=op.gradient{q}(:,:,ix); grad=pagemtimes(G,ue);
        strain=cat(1,grad(1,1,:,:),grad(2,2,:,:),grad(3,3,:,:), ...
            grad(1,2,:,:)+grad(2,1,:,:),grad(2,3,:,:)+grad(3,2,:,:),grad(1,3,:,:)+grad(3,1,:,:));
        sigma=pagemtimes(D,strain); S=zeros(3,3,ne,nm);
        S(1,1,:,:)=sigma(1,1,:,:); S(2,2,:,:)=sigma(2,1,:,:); S(3,3,:,:)=sigma(3,1,:,:);
        S(1,2,:,:)=sigma(4,1,:,:);S(2,1,:,:)=sigma(4,1,:,:);
        S(2,3,:,:)=sigma(5,1,:,:);S(3,2,:,:)=sigma(5,1,:,:);
        S(1,3,:,:)=sigma(6,1,:,:);S(3,1,:,:)=sigma(6,1,:,:);
        local=local+pagemtimes(permute(G,[2,1,3]),S).*reshape(op.weights(ix,q),1,1,ne,1);
    end
    values=permute(reshape(permute(local,[2,1,3,4]),3*na,ne,nm),[1,3,2]);
    values=values-pagemtimes(rigid,pagemtimes(permute(rigid,[2,1,3]),values));
    values=reshape(permute(values,[1,3,2]),[],nm);
    row=repmat(ids(:),nm,1);col=repelem((1:nm)',numel(ids));
    force=force+sparse(row,col,values(:),op.nd,nm);
end
end
