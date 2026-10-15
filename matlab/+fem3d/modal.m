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
free=bc.free; K=system.K(free,free); M=system.M(free,free);
assert(count<numel(free),'fem3d:Modes','Request fewer modes than free DOFs.');
[~,p,~]=chol(K); assert(p==0,'fem3d:Unconstrained','Modal supports must remove rigid body modes.');
if numel(free)<300
    [V,D]=eig(full(K),full(M)); [lambda,order]=sort(real(diag(D)));
    V=V(:,order(1:count)); lambda=lambda(1:count);
else
    opts=struct('tol',1e-9,'maxit',1500);
    [V,D,flag]=eigs(K,M,count,'smallestabs',opts);
    assert(flag==0,'fem3d:Modes','Eigenvalue iteration did not converge.');
    [lambda,order]=sort(real(diag(D))); V=real(V(:,order));
end
assert(all(lambda>0),'fem3d:Modes','Nonpositive structural eigenvalues.');
residual=zeros(count,1);
for k=1:count
    V(:,k)=V(:,k)/sqrt(V(:,k)'*M*V(:,k));
    [~,idx]=max(abs(V(:,k))); V(:,k)=V(:,k)*sign(V(idx,k));
    residual(k)=norm(K*V(:,k)-lambda(k)*M*V(:,k))/ ...
        max(norm(K*V(:,k))+norm(lambda(k)*M*V(:,k)),eps);
end
Phi=zeros(size(system.K,1),count); Phi(free,:)=V;
result=struct('frequencyHz',sqrt(lambda)/(2*pi),'omega',sqrt(lambda), ...
    'Phi',Phi,'eigenResidual',residual,'massNormalization',"Phi'*M*Phi=I",'bc',bc);
end
