function result = static(K,fixed,F,O)
%STATIC Full linear equilibrium with homogeneous supports and explicit outputs.
nd=size(K,1); bc=fem3d.prepare_constraints(nd,struct('dofs',fixed));
assert(size(F,1)==nd && all(isfinite(F(:))) && size(O,2)==nd, ...
    'ael:Load','Invalid force/output dimensions.');
f=bc.free; [L,p,order]=chol(K(f,f),'lower','vector');
assert(p==0,'ael:Constraints','Kff is not positive definite.');
U=zeros(nd,size(F,2)); U(f(order),:)=L'\(L\F(f(order),:)); R=K*U-F;
result=struct('U',U,'Y',O*U,'reaction',R(bc.dofs,:), ...
    'relativeResidual',norm(R(f,:),'fro')/max(norm(F(f,:),'fro'),eps));
end
