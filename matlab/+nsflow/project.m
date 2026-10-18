function [U,p,defect] = project(g,U,dt,rho)
%PROJECT Exact discrete MAC pressure projection by periodic FFT Poisson solve.
% The zero pressure mode is fixed to zero; no mass source is removed silently.
div=nsflow.divergence(g,U);rhs=rho/dt*div;
assert(abs(mean(rhs,'all'))<1e-9*max(1,max(abs(rhs),[],'all')),'nsflow:Mass','Periodic flux imbalance.');
hat=fftn(rhs)./g.laplacianEigenvalues;hat(1)=0;p=real(ifftn(hat));
for i=1:g.dimension,U{i}=U{i}-dt/rho*(p-nsflow.shift(p,i,1))/g.h(i);end
defect=max(abs(nsflow.divergence(g,U)),[],'all');
end
