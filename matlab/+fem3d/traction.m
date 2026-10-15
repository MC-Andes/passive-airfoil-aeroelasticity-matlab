function F = traction(mesh,faces,load)
%TRACTION Consistent nodal forces for constant traction/pressure per triangle.
% Scalar positive pressure acts INWARD on outward-oriented boundary faces.
% load: scalar pressure Pa; Nface x 1 pressures; 1 x 3 traction Pa;
% Nface x 3 tractions Pa. For a single face, vector 1x3 means traction.
nn=size(mesh.nodes,1); nf=size(faces,1); X=mesh.nodes;
normal=cross(X(faces(:,2),:)-X(faces(:,1),:),X(faces(:,3),:)-X(faces(:,1),:),2);
area=.5*vecnorm(normal,2,2);
assert(all(area>0) && all(isfinite(load(:))),'fem3d:Load','Invalid faces or load.');
if isscalar(load)
    total=-.5*load*normal;
elseif isequal(size(load),[1,3]) || isequal(size(load),[nf,3])
    total=area.*load;
else
    assert(isequal(size(load),[nf,1]),'fem3d:Load','Invalid pressure/traction dimensions.');
    total=-.5*load.*normal;
end
F=zeros(nn,3);
for component=1:3
    F(:,component)=accumarray(faces(:),repmat(total(:,component)/3,3,1),[nn,1]);
end
end
