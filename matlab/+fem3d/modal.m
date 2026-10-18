function result = modal(system,bc,count)
%MODAL Structural eigenmodes from assembled K/M; no aeroelastic flutter claim.
arguments
    system struct
    bc struct
    count (1,1) double {mustBeInteger,mustBePositive} = 6
end
bc=fem3d.prepare_constraints(size(system.K,1),bc);
fem3d.check_supports(system.mesh,bc);
assert(all(bc.values==0),'fem3d:BC','Modal analysis requires homogeneous BC.');
free=bc.free; Kphysical=system.K(free,free); Mphysical=system.M(free,free);
scale=1./sqrt(diag(Kphysical)); S=spdiags(scale,0,numel(free),numel(free));
K=S*Kphysical*S; M=S*Mphysical*S;
assert(count<numel(free),'fem3d:Modes','Request fewer modes than free DOFs.');
[order,orderingInfo]=fem3d.sparse_order(K); [L,p]=chol(K(order,order),'lower');
assert(p==0,'fem3d:Unconstrained','Modal supports must remove rigid body modes.');
if numel(free)<300
    [V,D]=eig(full(K),full(M)); [lambda,eigenOrder]=sort(real(diag(D)));
    V=V(:,eigenOrder(1:count)); lambda=lambda(1:count);
else
    opts=struct('tol',1e-10,'maxit',1500,'issym',true,'isreal',true);
    inverse=@(x)fem3d.factor_solve(L,order,x);
    [V,D,flag]=eigs(inverse,numel(free),M,count,0,opts);
    assert(flag==0,'fem3d:Modes','Eigenvalue iteration did not converge.');
    [lambda,eigenOrder]=sort(real(diag(D))); V=real(V(:,eigenOrder));
end
assert(all(lambda>0),'fem3d:Modes','Nonpositive structural eigenvalues.');
V=scale.*V; K=Kphysical; M=Mphysical;
residual=zeros(count,1);
for k=1:count
    V(:,k)=V(:,k)/sqrt(V(:,k)'*M*V(:,k));
    [~,idx]=max(abs(V(:,k))); V(:,k)=V(:,k)*sign(V(idx,k));
    residual(k)=norm(K*V(:,k)-lambda(k)*M*V(:,k))/ ...
        max(norm(K*V(:,k))+norm(lambda(k)*M*V(:,k)),eps);
end
assembledResidual=residual; residualDefinition='assembled sparse matrix'; stableHistory=[];
if isfield(system.mesh,'elevation') && max(residual)>1e-7
    op=fem3d.action_operator(system); fullV=zeros(size(system.K,1),count); fullV(free,:)=V;
    action=fem3d.action(op,fullV); action=action(free,:); oldMetric=Inf;bestV=V;bestAction=action;bestLambda=lambda;
    for iteration=1:10
        mv=M*V; defect=action-mv.*lambda';
        current=vecnorm(defect)./max(vecnorm(action)+vecnorm(mv.*lambda'),realmin);
        stableHistory(end+1,1)=max(current); %#ok<AGROW>
        if max(current)>=oldMetric,V=bestV;action=bestAction;lambda=bestLambda;break;end
        bestV=V;bestAction=action;bestLambda=lambda;
        if max(current)<1e-7,break;end
        oldMetric=max(current);
        correction=scale.*fem3d.factor_solve(L,order,scale.*defect);
        candidate=V-correction; R=chol(candidate'*M*candidate); candidate=candidate/R;
        fullV(free,:)=candidate; candidateAction=fem3d.action(op,fullV); candidateAction=candidateAction(free,:);
        reduced=candidate'*candidateAction; reduced=(reduced+reduced')/2;
        [rotation,eigenvalues]=eig(reduced); [lambdaNew,sortIndex]=sort(real(diag(eigenvalues)));rotation=rotation(:,sortIndex);
        V=candidate*rotation;action=candidateAction*rotation;lambda=lambdaNew;
    end
    for mode=1:count
        residual(mode)=norm(action(:,mode)-lambda(mode)*M*V(:,mode))/ ...
            max(norm(action(:,mode))+norm(lambda(mode)*M*V(:,mode)),realmin);
    end
    residualDefinition='elementwise strain/stress action with local rigid translation removed';
end
if ~exist('action','var'),action=K*V;end
defect=action-(M*V).*lambda'; correction=scale.*fem3d.factor_solve(L,order,scale.*defect);
energyResidual=sqrt(max(sum(defect.*correction,1),0)./max(sum(V.*action,1),realmin));
assembledResidual=zeros(count,1);
for mode=1:count
    assembledResidual(mode)=norm(K*V(:,mode)-lambda(mode)*M*V(:,mode))/ ...
        max(norm(K*V(:,mode))+norm(lambda(mode)*M*V(:,mode)),realmin);
end
Phi=zeros(size(system.K,1),count); Phi(free,:)=V;
result=struct('frequencyHz',sqrt(lambda)/(2*pi),'omega',sqrt(lambda), ...
    'Phi',Phi,'eigenResidual',residual,'eigenEnergyResidual',energyResidual(:),'assembledEigenResidual',assembledResidual, ...
    'residualDefinition',residualDefinition,'stableActionResidualHistory',stableHistory,'orderingInfo',orderingInfo,'factorNonzeros',nnz(L),'massNormalization',"Phi'*M*Phi=I",'bc',bc);
end
