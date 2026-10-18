function out = transform(axisOrder,origin)
%TRANSFORM Explicit CAD->solver convention for polar and axial vectors.
% Column vectors: x_solver=R*x_CAD-origin', moment_solver=det(R)*R*moment_CAD.
assert(isequal(sort(axisOrder),1:3),'cad:Axes','Expected an axis permutation.');
R=eye(3); R=R(axisOrder,:);
out=struct('R',R,'origin',origin,'determinant',det(R),'momentMap',det(R)*R, ...
    'positionRule','x_solver=R*x_CAD-origin; moment origin must be translated first', ...
    'materialBasisRule','transform constitutive tensor or rebuild a right-handed basis; no blind improper rotation');
end
