function F = traction(mesh,faces,load)
%TRACTION Consistent Tri3/Tri6 pressure or Cartesian traction.
% Positive scalar pressure acts INWARD. Numeric data: scalar, Nface x 1,
% 1 x 3 or Nface x 3. Callback @(points,normals) returns Nq x 1 or Nq x 3.
q=fem3d.surface_quadrature(mesh,faces); nf=size(faces,1); nq=numel(q.weights);
if isa(load,'function_handle')
    values=load(q.points,q.normals);
    assert(isequal(size(values),[nq,1]) || isequal(size(values),[nq,3]), ...
        'fem3d:Load','Callback must return Nq x 1 pressure or Nq x 3 traction.');
elseif isscalar(load)
    values=repmat(load,nq,1);
elseif isequal(size(load),[1,3])
    values=repmat(load,nq,1);
elseif isequal(size(load),[nf,1]) || isequal(size(load),[nf,3])
    values=load(q.faceIndex,:);
else
    error('fem3d:Load','Invalid pressure/traction dimensions.');
end
assert(all(isfinite(values(:))),'fem3d:Load','Nonfinite pressure/traction.');
if size(values,2)==1, values=-values.*q.normals; end
F=q.H'*(q.weights.*values);
end
