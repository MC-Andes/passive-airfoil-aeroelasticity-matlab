function [H,info] = transfer(nodes,points,candidates,neighbors)
%TRANSFER Local affine displacement interpolation and conservative force map.
% H maps structural nodal displacements N x 3 -> point displacements P x 3.
% Nodal forces = H' * point forces: exact force, moment and virtual work when
% affine reproduction succeeds. Not a CFD pressure reconstruction.
% User must select candidates belonging to the correct physical moving skin.
% Unrestricted nearest-neighbor patches may span a gap or unrelated component.
arguments
    nodes (:,3) double {mustBeFinite}
    points (:,3) double {mustBeFinite}
    candidates double {mustBeInteger,mustBePositive} = (1:size(nodes,1))'
    neighbors (1,1) double {mustBeInteger,mustBePositive} = 16
end
candidates=unique(candidates(:)); nn=size(nodes,1); np=size(points,1);
assert(numel(candidates)>=4 && all(candidates<=nn),'fem3d:Transfer','Insufficient/invalid transfer nodes.');
I=cell(np,1); J=I; W=I; err=zeros(np,1); radius=err; negative=err;
for p=1:np
    dist2=sum((nodes(candidates,:)-points(p,:)).^2,2);
    [~,order]=sort(dist2); k=min(neighbors,numel(order));
    while true
        ids=candidates(order(1:k)); relative=nodes(ids,:)-points(p,:);
        scale=max(vecnorm(relative,2,2)); scale=max(scale,eps);
        A=[ones(1,k);relative'/scale]; b=[1;0;0;0];
        % Minimum weighted norm subject to affine reproduction; SVD also
        % handles a planar patch when the target lies in that plane.
        sw=1./sqrt(1+sum(relative.^2,2)/scale^2);
        [U,S,V]=svd(A.*sw','econ'); s=diag(S);
        keep=s>1e-11*s(1);
        w=sw.*(V(:,keep)*((U(:,keep)'*b)./s(keep)));
        defect=norm(A*w-b,inf);
        if defect<1e-9, break; end
        if k==numel(order) || k>=256
            error('fem3d:Transfer','Cannot reproduce affine motion at load point %d; review skin patch.',p);
        end
        k=min([2*k,numel(order),256]);
    end
    I{p}=repmat(p,k,1); J{p}=ids; W{p}=w;
    err(p)=defect; radius(p)=max(vecnorm(relative,2,2)); negative(p)=sum(w<0);
end
H=sparse(vertcat(I{:}),vertcat(J{:}),vertcat(W{:}),np,nn);
info=struct('maximumAffineDefect',max(err),'maximumPatchRadius',max(radius), ...
    'negativeWeightCount',sum(negative),'patchRadii',radius, ...
    'method','local affine minimum-norm interpolation; force transpose by virtual work');
end
