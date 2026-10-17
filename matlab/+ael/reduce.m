function model = reduce(K,M,fixed,count,zeta,B,O)
%REDUCE Linear FEM modal model with explicit load and measurement operators.
% Works with solid or plate matrices; never assumes their DOF conventions.
% B maps input amplitudes to generalized nodal forces; O maps DOFs to outputs.
% Homogeneous supports, mass-normalized modes, supplied modal damping ratios.
nd=size(K,1); fixed=unique(fixed(:)); free=setdiff((1:nd)',fixed);
assert(isequal(size(K),size(M),[nd,nd]) && all(isfinite(nonzeros(K))) && ...
    all(isfinite(nonzeros(M))) && norm(K-K','fro')<1e-10*max(norm(K,'fro'),eps) && ...
    norm(M-M','fro')<1e-10*max(norm(M,'fro'),eps),'ael:Matrices','Invalid symmetric K/M.');
assert(all(fixed>=1 & fixed<=nd & fixed==fix(fixed)) && count>=1 && ...
    count==fix(count) && count<numel(free),'ael:Constraints','Invalid supports or mode count.');
assert(size(B,1)==nd && size(O,2)==nd && all(isfinite(nonzeros(B))) && ...
    all(isfinite(nonzeros(O))),'ael:Operators','Incompatible input/output operators.');
if isscalar(zeta), zeta=repmat(zeta,count,1); end
assert(numel(zeta)==count && all(isfinite(zeta(:)) & zeta(:)>=0), ...
    'ael:Damping','Supply nonnegative damping per retained mode.');
Kf=K(free,free); Mf=M(free,free);
[~,pk,~]=chol(Kf); [~,pm,~]=chol(Mf);
assert(pk==0 && pm==0,'ael:Constraints','Kff/Mff must be positive definite.');
if numel(free)<250
    [V,D]=eig(full(Kf),full(Mf)); [lambda,order]=sort(real(diag(D)));
    V=V(:,order(1:count)); lambda=lambda(1:count);
else
    [V,D,flag]=eigs(Kf,Mf,count,'smallestabs',struct('tol',1e-10,'maxit',2000));
    assert(flag==0,'ael:Modes','Eigenvalue iteration did not converge.');
    [lambda,order]=sort(real(diag(D))); V=real(V(:,order));
end
assert(all(lambda>0),'ael:Modes','Nonpositive eigenvalues.');
residual=zeros(count,1);
for k=1:count
    V(:,k)=V(:,k)/sqrt(V(:,k)'*Mf*V(:,k));
    [~,idx]=max(abs(V(:,k))); V(:,k)=V(:,k)*sign(V(idx,k));
    residual(k)=norm(Kf*V(:,k)-lambda(k)*Mf*V(:,k))/ ...
        (norm(Kf*V(:,k))+norm(lambda(k)*Mf*V(:,k)));
end
Phi=zeros(nd,count); Phi(free,:)=V; omega=sqrt(lambda);
model=struct('M',Phi'*M*Phi,'K',Phi'*K*Phi, ...
    'C',diag(2*zeta(:).*omega),'B',Phi'*B,'O',O*Phi,'Phi',Phi, ...
    'frequencyHz',omega/(2*pi),'omega',omega,'zeta',zeta(:),'free',free, ...
    'eigenResidual',residual,'scope','decoupled linear structural modal model');
end
