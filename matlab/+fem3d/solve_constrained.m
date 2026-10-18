function r = solve_constrained(system,F,bc,C,d)
%SOLVE_CONSTRAINED Exact affine multipoint constraints C*u=d via multipliers.
% Supports are eliminated first. Disconnected physical parts may be joined by
% documented C rows; no stabilization springs and no inferred CAD interfaces.
nd=size(system.K,1); if size(F,2)==3,F=reshape(F',[],1);else,F=F(:);end
bc=fem3d.prepare_constraints(nd,bc);free=bc.free;fixed=bc.dofs;
assert(size(C,2)==nd && numel(d)==size(C,1),'fem3d:MPC','Invalid constraint dimensions.');
Cf=C(:,free);d=d(:)-C(:,fixed)*bc.values; rhs=F(free)-system.K(free,fixed)*bc.values;
K=system.K(free,free); nc=size(C,1); scale=1./sqrt(max(diag(K),realmin));S=spdiags(scale,0,numel(free),numel(free));
rowScale=1./sqrt(sum((Cf*S).^2,2));assert(all(isfinite(rowScale)),'fem3d:MPC','Redundant constraints on prescribed DOFs.');
R=spdiags(rowScale,0,nc,nc); A=[S*K*S,S*Cf'*R;R*Cf*S,sparse(nc,nc)];
warningState=warning('off','MATLAB:nearlySingularMatrix');cleanup=onCleanup(@()warning(warningState));
x=A\[scale.*rhs;rowScale.*d];u=zeros(nd,1);u(fixed)=bc.values;u(free)=scale.*x(1:numel(free));lambda=rowScale.*x(numel(free)+1:end);
constraint=C*u-(d+C(:,fixed)*bc.values);residual=system.K*u+C'*lambda-F;
backward=norm(residual(free),inf)/max(norm(K,1)*norm(u(free),inf)+norm(F,inf),realmin);
assert(all(isfinite(x)) && (isempty(constraint) || norm(constraint,inf)<1e-9) && backward<1e-10, ...
    'fem3d:MPC','Singular/inconsistent multipoint system.');
r=struct('u',u,'displacement',reshape(u,3,[])','lambda',lambda,'constraintResidual',constraint, ...
    'equilibriumResidual',residual,'normwiseBackwardError',backward,'reactionFromMPC',-C'*lambda,'bc',bc);
end
