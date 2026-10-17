function result = newmark_observed(system,bc,time,forceFcn,O,u0,v0,rayleigh)
%NEWMARK_OBSERVED Full-order Newmark with streamed states and explicit sensors.
% Avoids storing 3N x Nt displacement/velocity/acceleration histories. Matrices
% and current full states are retained; no modal reduction is performed here.
if nargin<8,rayleigh=[0,0];end
nd=size(system.K,1);bc=fem3d.prepare_constraints(nd,bc);assert(all(bc.values==0),'fem3d:BC','Streamed comparison uses homogeneous supports.');
f=bc.free;time=time(:)';nt=numel(time);assert(nt>=2 && all(diff(time)>0) && size(O,2)==nd,'fem3d:Time','Invalid full-order observation/time.');
K=system.K(f,f);M=system.M(f,f);C=rayleigh(1)*M+rayleigh(2)*K;u=u0(f);v=v0(f);
load=forceFcn(time(1));if size(load,2)==3,load=reshape(load',[],1);else,load=load(:);end
rhs=load(f)-C*v-K*u;[a,flag]=pcg(M,rhs,1e-11,500,spdiags(diag(M),0,numel(f),numel(f)));assert(flag==0,'fem3d:Mass','Initial consistent-mass iteration failed.');
Y=zeros(size(O,1),nt);energy=zeros(1,nt);work=energy;diss=energy;residual=energy;Y(:,1)=O*u0;
energy(1)=.5*(v'*M*v+u'*K*u);previousLoad=load(f);oldDt=NaN;
for i=2:nt
    dt=time(i)-time(i-1);
    if ~isfinite(oldDt) || abs(dt-oldDt)>1e-10*dt
        effective=M+.5*dt*C+.25*dt^2*K;scale=1./sqrt(diag(effective));S=spdiags(scale,0,numel(f),numel(f));A=S*effective*S;
        order=fem3d.sparse_order(A);[L,p]=chol(A(order,order),'lower');assert(p==0,'fem3d:Time','Full effective matrix is not SPD.');oldDt=dt;
    end
    up=u+dt*v+.25*dt^2*a;vp=v+.5*dt*a;
    load=forceFcn(time(i));if size(load,2)==3,load=reshape(load',[],1);else,load=load(:);end
    rhs=load(f)-C*vp-K*up;nextA=scale.*fem3d.factor_solve(L,order,scale.*rhs);
    nextU=up+.25*dt^2*nextA;nextV=vp+.5*dt*nextA;
    Y(:,i)=O(:,f)*nextU;energy(i)=.5*(nextV'*M*nextV+nextU'*K*nextU);
    du=nextU-u;vmid=(nextV+v)/2;work(i)=work(i-1)+.5*(load(f)+previousLoad)'*du;diss(i)=diss(i-1)+dt*vmid'*C*vmid;
    residual(i)=norm(M*nextA+C*nextV+K*nextU-load(f),inf);
    u=nextU;v=nextV;a=nextA;previousLoad=load(f);
end
balance=energy-energy(1)-work+diss;
result=struct('time',time,'Y',Y,'energy',energy,'externalWork',work,'dissipation',diss,'energyBalance',balance, ...
    'relativeEnergyBalance',max(abs(balance))/max([max(abs(energy)),max(abs(work)),1e-15]), ...
    'maximumForceResidual',max(residual),'lastU',u,'lastV',v,'rayleigh',rayleigh,'factorNonzeros',nnz(L), ...
    'scope','streamed FULL solid FEM Newmark, no modal truncation');
end
