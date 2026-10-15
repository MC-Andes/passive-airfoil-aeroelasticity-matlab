function [nodes,J] = geometry(model,q)
%GEOMETRY Geometry and virtual-displacement Jacobian at current coordinates.
q=q(:); assert(numel(q)==model.ndof && all(isfinite(q)),'pae:State','Invalid coordinates.');
X0=reshape(model.surface.nodes,[],3);
if isfield(model.structure,'geometryFcn')
    [X,J]=model.structure.geometryFcn(q);
else
    J=model.structure.Phi;
    X=reshape(X0(:)+J*q,size(X0));
end
assert(isequal(size(X),size(X0)) && isequal(size(J),[numel(X0),model.ndof]) ...
    && all(isfinite([X(:);J(:)])),'pae:Geometry','Invalid geometryFcn coordinates/Jacobian.');
nodes=reshape(X,size(model.surface.nodes));
end
