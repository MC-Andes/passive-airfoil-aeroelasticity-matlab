function Phi = interpolate_basis(samplePoints,sampleModes,targetPoints)
%INTERPOLATE_BASIS Cubic RBF with affine reproduction for sampled ROM modes.
% sampleModes: N-by-3-by-r; Phi: (3*Ntarget)-by-r in MATLAB X(:) ordering.
% Use a modest set of surface samples; solves a dense interpolation system.
assert(size(samplePoints,2)==3 && size(targetPoints,2)==3 && ...
    size(sampleModes,1)==size(samplePoints,1) && size(sampleModes,2)==3 && ...
    all(isfinite([samplePoints(:);sampleModes(:);targetPoints(:)])), ...
    'pae:Interpolation','Invalid modal sample data.');
n=size(samplePoints,1);
assert(n>=3 && n<=1500,'pae:Interpolation','Use 3 to 1500 representative RBF samples.');
assert(size(unique(samplePoints,'rows'),1)==n,'pae:Interpolation','Duplicate modal sample points.');
origin=mean(samplePoints,1); centered=samplePoints-origin;
[~,S,V]=svd(centered,0); singular=diag(S);
dimension=sum(singular>max(singular)*1e-10);
assert(dimension>=1,'pae:Interpolation','Degenerate modal sample coordinates.');
scale=max(vecnorm(centered,2,2));
source=centered*V(:,1:dimension)/scale;
target=(targetPoints-origin)*V(:,1:dimension)/scale;
outside=(targetPoints-origin)-target*scale*V(:,1:dimension)';
assert(max(vecnorm(outside,2,2))<1e-7*max(1,scale), ...
    'pae:Interpolation','Targets lie outside the sample affine subspace.');
P=[ones(n,1),source]; p=dimension+1;
D=distance(source,source).^3;
system=[D,P;P',zeros(p)];
assert(rcond(system)>1e-14,'pae:Interpolation','Ill-conditioned RBF samples.');
r=size(sampleModes,3); values=reshape(sampleModes,n,3*r);
coefficients=system\[values;zeros(p,3*r)];
H=[distance(target,source).^3,ones(size(target,1),1),target];
interpolated=reshape(H*coefficients,size(target,1),3,r);
Phi=reshape(interpolated,3*size(target,1),r);
end

function D=distance(a,b)
D=sqrt(max(0,sum(a.^2,2)+sum(b.^2,2)'-2*a*b'));
end
