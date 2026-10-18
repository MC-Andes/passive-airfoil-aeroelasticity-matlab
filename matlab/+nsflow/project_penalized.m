function [U,p,info] = project_penalized(g,R,t,dt,physical)
%PROJECT_PENALIZED Coupled implicit Brinkman and variable-coefficient pressure.
% Unlike split penalty/constant projection, the pressure correction includes
% the wall coefficient, so projection does not reintroduce a dt-sized wall
% slip. PCG on the mean-zero space uses the own periodic FFT preconditioner.
d=g.dimension;beta=cell(1,d);wall=beta;mask=beta;U=beta;
for a=1:d
 mask{a}=physical.mask(g.velocityCoordinates(a,:),t);wall{a}=zeros(g.count);if isfield(physical,'wallVelocity'),wall{a}=physical.wallVelocity(g.velocityCoordinates(a,:),t,a);end
 beta{a}=1./(1+dt*mask{a}/physical.penaltyTime);U{a}=beta{a}.*(R{a}+dt*mask{a}.*wall{a}/physical.penaltyTime);
end
rhs=-physical.rho/dt*nsflow.divergence(g,U);rhs=rhs(:);assert(abs(mean(rhs))<1e-9*max(norm(rhs,inf),1),'nsflow:Mass','Nonzero periodic source.');
if norm(rhs)==0,p=zeros(g.count);flag=0;relres=0;iteration=0;
else,[p,flag,relres,iteration]=pcg(@apply,rhs,1e-10,150,@precondition);p=reshape(p,g.count);p=p-mean(p,'all');assert(flag==0,'nsflow:Pressure','Variable-coefficient pressure iteration failed: %d, residual %g.',flag,relres);end
body=zeros(1,d);work=0;
for a=1:d
 U{a}=U{a}-dt/physical.rho*beta{a}.*(p-nsflow.shift(p,a,1))/g.h(a);force=physical.rho*mask{a}.*(U{a}-wall{a})/physical.penaltyTime;body(a)=g.cellVolume*sum(force,'all');work=work-g.cellVolume*sum(force.*wall{a},'all');
end
info=struct('bodyForce',body,'wallWork',work,'maximumDivergence',max(abs(nsflow.divergence(g,U)),[],'all'),'pressureIterations',iteration,'pressureRelativeResidual',relres);
 function y=apply(x)
  x=reshape(x,g.count);gradient=cell(1,d);for j=1:d,gradient{j}=beta{j}.*(x-nsflow.shift(x,j,1))/g.h(j);end;y=-nsflow.divergence(g,gradient);y=y(:);
 end
 function y=precondition(x)
  x=reshape(x,g.count);hat=-fftn(x)./g.laplacianEigenvalues;hat(1)=0;y=real(ifftn(hat));y=y(:);
 end
end
