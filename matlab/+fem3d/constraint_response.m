function response = constraint_response(system,F,bc,C)
%CONSTRAINT_RESPONSE One factorization for an exact guide compliance/Schur map.
% Returns Kff^-1 F and Kff^-1 C' with energy-residual estimates. This requires
% an SPD monolithic base; solve_constrained retains KKT for other mechanisms.
nd=size(system.K,1);if isequal(size(F),[nd/3,3]),F=reshape(F',[],1);end
assert(size(F,1)==nd,'fem3d:Guide','Loads must be 3N by load count or N by 3.');
bc=fem3d.prepare_constraints(nd,bc);assert(all(bc.values==0),'fem3d:Guide','This family requires homogeneous base supports.');
f=bc.free;K=system.K(f,f);scale=full(1./sqrt(diag(K)));S=spdiags(scale,0,numel(f),numel(f));A=S*K*S;
[order,info]=fem3d.sparse_order(A);[L,status]=chol(A(order,order),'lower');assert(status==0,'fem3d:Guide','Base structure must be SPD.');
rhs=full([F(f,:),full(C(:,f)')]);U=zeros(nd,size(rhs,2));U(f,:)=scale.*fem3d.factor_solve(L,order,scale.*rhs);
op=fem3d.action_operator(system);energyResidual=zeros(1,size(U,2));
for iteration=1:3
 internal=fem3d.action(op,U);r=rhs-internal(f,:);correction=scale.*fem3d.factor_solve(L,order,scale.*r);
 energyResidual=sqrt(max(sum(r.*correction,1),0)./max(sum(U.*internal,1),realmin));
 if max(energyResidual)<1e-8,break;end
 U(f,:)=U(f,:)+correction;
end
nf=size(F,2);response=struct('base',U(:,1:nf),'influence',U(:,nf+1:end),'compliance',C*U(:,nf+1:end), ...
 'observation',C,'energyResidual',energyResidual,'bc',bc,'orderingInfo',info);
response.compliance=(response.compliance+response.compliance')/2;
end
