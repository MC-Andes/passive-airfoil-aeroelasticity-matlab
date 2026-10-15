function bc = prepare_constraints(nd,bc)
%PREPARE_CONSTRAINTS Exact elimination of constant prescribed displacements.
assert(isfield(bc,'dofs') && ~isempty(bc.dofs),'fem3d:BC','Specify constrained DOFs.');
fixed=bc.dofs(:);
assert(all(fixed==fix(fixed) & fixed>=1 & fixed<=nd) && numel(unique(fixed))==numel(fixed), ...
    'fem3d:BC','Constrained DOFs must be unique and valid.');
values=zeros(numel(fixed),1); if isfield(bc,'values'), values=bc.values(:); end
assert(numel(values)==numel(fixed) && all(isfinite(values)),'fem3d:BC','Invalid BC values.');
bc.dofs=fixed; bc.values=values; bc.free=setdiff((1:nd)',fixed);
assert(~isempty(bc.free),'fem3d:BC','No free DOFs.');
end
