function [force,tangent] = internal_force(model,q)
%INTERNAL_FORCE Custom constitutive law or linear plus diagonal cubic law.
q=q(:);
if isfield(model.structure,'internalFcn')
    [force,tangent]=model.structure.internalFcn(q);
else
    force=model.structure.K*q+model.structure.cubic.*q.^3;
    tangent=model.structure.K+diag(3*model.structure.cubic.*q.^2);
end
assert(isequal(size(force),[model.ndof,1]) && isequal(size(tangent),[model.ndof,model.ndof]) ...
    && all(isfinite([force(:);tangent(:)])),'pae:Structure','Invalid restoring force/tangent.');
end
