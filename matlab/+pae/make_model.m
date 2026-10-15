function model = make_model(surface,structure)
%MAKE_MODEL Structural coordinates and aerodynamic surface are independent.
% Required structure: M,C,K (r-by-r), Phi (3*Nnodes-by-r), labels (r entries).
% X(:)=X0(:)+Phi*q uses MATLAB column order: all x, then all y, then all z.
% Optional: cubic (r entries), externalForce, internalFcn(q)->[force,tangent],
% geometryFcn(q)->[N-by-3 coordinates,3*N-by-r Jacobian], aeroFcn(nodes,V,flow).
% The builtin structure is a constant-mass reduced model, not nonlinear FEM.
assert(isfield(surface,'nodes'),'pae:Model','Missing surface nodes.');
pae.panels(surface.nodes);
required={'M','C','K','Phi','labels','source'};
for i=1:numel(required)
    assert(isfield(structure,required{i}),'pae:Model','Missing structural field %s.',required{i});
end
r=size(structure.M,1); N=numel(surface.nodes)/3;
assert(r>=1 && isequal(size(structure.M),[r,r]) && isequal(size(structure.C),[r,r]) ...
    && isequal(size(structure.K),[r,r]) && isequal(size(structure.Phi),[3*N,r]) ...
    && numel(structure.labels)==r,'pae:Model','Incompatible structural dimensions.');
for name={'M','C','K'}
    matrix=structure.(name{1});
    assert(all(isfinite(matrix(:))) && norm(matrix-matrix','fro')<1e-10*max(1,norm(matrix,'fro')), ...
        'pae:Model','%s must be finite and symmetric.',name{1});
end
[~,flag]=chol(structure.M); assert(flag==0,'pae:Model','Mass must be positive definite.');
assert(min(eig(structure.K))>=-1e-10*max(1,norm(structure.K)), ...
    'pae:Model','Builtin stiffness must be positive semidefinite. Use internalFcn for other laws.');
assert(min(eig(structure.C))>=-1e-10*max(1,norm(structure.C)), ...
    'pae:Model','Passive structural damping must be positive semidefinite.');
assert(all(isfinite(structure.Phi(:))),'pae:Model','Invalid deformation basis.');
if ~isfield(structure,'cubic'), structure.cubic=zeros(r,1); end
if ~isfield(structure,'externalForce'), structure.externalForce=zeros(r,1); end
structure.cubic=structure.cubic(:); structure.externalForce=structure.externalForce(:);
if ~isfield(structure,'coordinateUnits')
    structure.coordinateUnits=repmat({'unidades del usuario'},r,1);
end
assert(numel(structure.coordinateUnits)==r && all(strlength(string(structure.coordinateUnits))>0), ...
    'pae:Model','Specify one coordinate unit per structural coordinate.');
structure.coordinateUnits=cellstr(string(structure.coordinateUnits(:)));
assert(numel(structure.cubic)==r && numel(structure.externalForce)==r && ...
    all(isfinite([structure.cubic;structure.externalForce])), ...
    'pae:Model','Invalid cubic/external force data.');
model=struct('surface',surface,'structure',structure,'ndof',r,'aeroOptions',struct());
end
