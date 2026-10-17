function result = solve_static(system,F,bc)
%SOLVE_STATIC Linear FEM equilibrium; constrained DOFs eliminated exactly.
% F is N x 3 nodal forces [N] or an interleaved 3N x 1 vector.
nn=size(system.mesh.nodes,1); nd=3*nn;
if isequal(size(F),[nn,3]), F=reshape(F',[],1); else, F=F(:); end
assert(numel(F)==nd && all(isfinite(F)),'fem3d:Load','Invalid nodal force.');
bc=fem3d.prepare_constraints(nd,bc); free=bc.free; fixed=bc.dofs;
fem3d.check_supports(system.mesh,bc);
% Three-output sparse Cholesky enables AMD ordering and controls fill-in.
Kff=system.K(free,free);
% Symmetric diagonal equilibration preserves energy and prescribed supports.
scale=1./sqrt(diag(Kff)); S=spdiags(scale,0,numel(free),numel(free));
scaledK=S*Kff*S; [order,orderingInfo]=fem3d.sparse_order(scaledK);
[L,p]=chol(scaledK(order,order),'lower');
assert(p==0,'fem3d:Unconstrained', ...
    'Kff is not positive definite: inspect disconnected parts, supports or mesh quality.');
u=zeros(nd,1); u(fixed)=bc.values;
rhs=F(free)-system.K(free,fixed)*bc.values;
scaled=scale.*rhs; y=zeros(numel(free),1); y(order)=L'\(L\scaled(order));
x=scale.*y; residualHistory=zeros(8,1); iterations=0;
% Refinement uses the ORIGINAL operator. Retain the best actual residual;
% stop at the numerical floor rather than claiming a converged solve.
for iteration=1:8
    defect=rhs-Kff*x; residualHistory(iteration)=norm(defect,inf)/max(norm(rhs,inf),eps);
    iterations=iteration;
    if residualHistory(iteration)<1e-8 || iteration==8, break; end
    scaled=scale.*defect; y(order)=L'\(L\scaled(order)); candidate=x+scale.*y;
    if norm(rhs-Kff*candidate,inf)>=norm(defect,inf), break; end
    x=candidate;
end
u(free)=x;
assembledResidual=system.K*u-F; internal=system.K*u; stableHistory=[];
residualDefinition='assembled sparse matrix';
if isfield(system.mesh,'elevation') && norm(assembledResidual(free),inf)/max(norm(rhs,inf),eps)>1e-7
    op=fem3d.action_operator(system); internal=fem3d.action(op,u);
    residualDefinition='elementwise strain/stress action with local rigid translation removed';
    for iteration=1:10
        defect=F(free)-internal(free); stableHistory(end+1,1)=norm(defect,inf)/max(norm(rhs,inf),eps); %#ok<AGROW>
        if stableHistory(end)<1e-7,break;end
        correction=scale.*fem3d.factor_solve(L,order,scale.*defect);
        candidate=u; candidate(free)=candidate(free)+correction; nextInternal=fem3d.action(op,candidate);
        if norm(F(free)-nextInternal(free),inf)>=norm(defect,inf),break;end
        u=candidate;internal=nextInternal;
    end
end
defect=F(free)-internal(free); correction=scale.*fem3d.factor_solve(L,order,scale.*defect);
energyDefect=max(real(defect'*correction),0); energyScale=max(real(u'*internal),realmin);
energyResidual=sqrt(energyDefect/energyScale);
residual=internal-F; reaction=zeros(nd,1); reaction(fixed)=residual(fixed);
U=reshape(u,3,[])'; R=reshape(reaction,3,[])'; X=system.mesh.nodes;
balance=sum(F_as_matrix(F)+R,1);
moment=sum(cross(X,F_as_matrix(F)+R,2),1);
result=struct('u',u,'displacement',U,'F',F,'reaction',R,'bc',bc, ...
    'freeResidual',norm(residual(free),inf), ...
    'relativeResidual',norm(residual(free),inf)/max(norm(rhs,inf),eps), ...
    'forceBalance',balance,'momentBalance',moment, ...
    'strainEnergy',.5*u'*internal,'assembledRelativeResidual', ...
    norm((system.K(free,:)*u-F(free)),inf)/max(norm(rhs,inf),eps), ...
    'residualDefinition',residualDefinition,'energyResidual',energyResidual, ...
    'estimatedDisplacementCorrection',norm(correction,inf),'stableActionResidualHistory',stableHistory,'maximumDisplacement',max(vecnorm(U,2,2)), ...
    'factorNonzeros',nnz(L),'orderingInfo',orderingInfo,'equilibrated',true, ...
    'refinementResidualHistory',residualHistory(1:iterations),'stress',fem3d.stress(system,u));
end

function out=F_as_matrix(F)
out=reshape(F,3,[])';
end
