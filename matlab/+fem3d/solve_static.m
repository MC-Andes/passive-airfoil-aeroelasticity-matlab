function result = solve_static(system,F,bc)
%SOLVE_STATIC Linear FEM equilibrium; constrained DOFs eliminated exactly.
% F is N x 3 nodal forces [N] or an interleaved 3N x 1 vector.
nn=size(system.mesh.nodes,1); nd=3*nn;
if isequal(size(F),[nn,3]), F=reshape(F',[],1); else, F=F(:); end
assert(numel(F)==nd && all(isfinite(F)),'fem3d:Load','Invalid nodal force.');
bc=fem3d.prepare_constraints(nd,bc); free=bc.free; fixed=bc.dofs;
fem3d.check_supports(system.mesh,bc);
% Three-output sparse Cholesky enables AMD ordering and controls fill-in.
Kff=system.K(free,free); [L,p,order]=chol(Kff,'lower','vector');
assert(p==0,'fem3d:Unconstrained', ...
    'Kff is not positive definite: inspect disconnected parts, supports or mesh quality.');
u=zeros(nd,1); u(fixed)=bc.values;
rhs=F(free)-system.K(free,fixed)*bc.values;
u(free(order))=L'\(L\rhs(order));
residual=system.K*u-F; reaction=zeros(nd,1); reaction(fixed)=residual(fixed);
U=reshape(u,3,[])'; R=reshape(reaction,3,[])'; X=system.mesh.nodes;
balance=sum(F_as_matrix(F)+R,1);
moment=sum(cross(X,F_as_matrix(F)+R,2),1);
result=struct('u',u,'displacement',U,'F',F,'reaction',R,'bc',bc, ...
    'freeResidual',norm(residual(free),inf), ...
    'relativeResidual',norm(residual(free),inf)/max(norm(rhs,inf),eps), ...
    'forceBalance',balance,'momentBalance',moment, ...
    'strainEnergy',.5*u'*system.K*u,'maximumDisplacement',max(vecnorm(U,2,2)), ...
    'factorNonzeros',nnz(L),'stress',fem3d.stress(system,u));
end

function out=F_as_matrix(F)
out=reshape(F,3,[])';
end
