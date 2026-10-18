function result = modal_constrained(system,bc,C,count)
%MODAL_CONSTRAINED Full FEM eigenmodes with exact affine homogeneous guides.
% Projected stiffness inverse enforces C*u=0 through the compliance Schur
% complement. Constraints add reactions, not springs or stiffness tuning.
nd=size(system.K,1);bc=fem3d.prepare_constraints(nd,bc);assert(all(bc.values==0));f=bc.free;K=system.K(f,f);M=system.M(f,f);Cf=C(:,f);scale=full(1./sqrt(diag(K)));S=spdiags(scale,0,numel(f),numel(f));A=S*K*S;
[order,info]=fem3d.sparse_order(A);[L,status]=chol(A(order,order),'lower');assert(status==0,'fem3d:Guide','SPD monolithic base required.');
Y=base_inverse(full(Cf'));G=Cf*Y;G=(G+G')/2;assert(rcond(full(G))>1e-14,'fem3d:Guide','Dependent guide coordinates.');
opts=struct('tol',1e-9,'maxit',500,'issym',true,'isreal',true);[V,D,flag]=eigs(@inverse,numel(f),M,count,0,opts);assert(flag==0,'fem3d:Modes','Constrained eigensolver did not converge.');[lambda,p]=sort(real(diag(D)));V=real(V(:,p));
Phi=zeros(nd,count);energy=zeros(count,1);constraint=energy;forward=energy;
for mode=1:count
 v=V(:,mode);v=v/sqrt(v'*M*v);Phi(f,mode)=v;r=K*v-lambda(mode)*M*v;reaction=(Cf*Cf')\(Cf*r);r=r-Cf'*reaction;
 correction=inverse(r);energy(mode)=sqrt(max(real(r'*correction),0)/max(real(v'*K*v),realmin));constraint(mode)=norm(Cf*v,inf);forward(mode)=norm(r)/max(norm(K*v)+norm(lambda(mode)*M*v),realmin);
end
result=struct('frequencyHz',sqrt(lambda)/(2*pi),'Phi',Phi,'eigenEnergyResidual',energy,'constraintResidual',constraint,'forwardResidual',forward,'orderingInfo',info,'scope','full monolithic FEM with exact mean-coordinate affine guide; conditional mechanism hypothesis');
 function x=base_inverse(b),x=scale.*fem3d.factor_solve(L,order,scale.*b);end
 function x=inverse(b),x=base_inverse(b);x=x-Y*(G\(Cf*x));end
end
