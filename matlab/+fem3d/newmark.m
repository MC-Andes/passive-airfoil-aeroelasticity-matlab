function result = newmark(system,bc,time,forceFcn,u0,v0,rayleigh)
%NEWMARK Own average-acceleration integration of assembled structural system.
% Decoupled externally specified nodal forces. Constant BC only.
% rayleigh=[mass coefficient 1/s, stiffness coefficient s], default zero.
arguments
    system struct
    bc struct
    time double
    forceFcn (1,1) function_handle
    u0 double = []
    v0 double = []
    rayleigh (1,2) double {mustBeNonnegative,mustBeFinite} = [0,0]
end
time=time(:)'; nt=numel(time); nd=size(system.K,1);
assert(nt>=2 && all(isfinite(time)) && all(diff(time)>0),'fem3d:Time','Require increasing times.');
bc=fem3d.prepare_constraints(nd,bc); f=bc.free; c=bc.dofs;
fem3d.check_supports(system.mesh,bc);
U=zeros(nd,nt); V=U; A=U;
if isempty(u0), u0=zeros(nd,1); u0(c)=bc.values; end
if isempty(v0), v0=zeros(nd,1); end
u0=u0(:); v0=v0(:);
assert(numel(u0)==nd && numel(v0)==nd && all(isfinite([u0;v0])) && ...
    norm(u0(c)-bc.values,inf)<1e-12 && norm(v0(c),inf)<1e-12, ...
    'fem3d:Initial','Initial state must satisfy constant BC.');
K=system.K(f,f); M=system.M(f,f); C=rayleigh(1)*M+rayleigh(2)*K;
[LM,p,morder]=chol(M,'lower','vector'); assert(p==0,'fem3d:Mass','Mass must be positive definite.');
[~,p,~]=chol(K); assert(p==0,'fem3d:Unconstrained','Insufficient supports.');
U(:,1)=u0; V(:,1)=v0;
force=getforce(forceFcn,time(1),nd)-system.K(:,c)*bc.values;
rhs=force(f)-C*v0(f)-K*u0(f);
A(f(morder),1)=LM'\(LM\rhs(morder));
beta=.25; gamma=.5; olddt=NaN; energy=zeros(1,nt); residual=zeros(1,nt);
energy(1)=.5*v0(f)'*M*v0(f)+.5*u0'*system.K*u0;
for i=2:nt
    dt=time(i)-time(i-1);
    if ~isfinite(olddt) || abs(dt-olddt)>1e-12*dt
        [L,p,order]=chol(M+gamma*dt*C+beta*dt^2*K,'lower','vector');
        assert(p==0,'fem3d:Time','Effective matrix is not positive definite.');
        olddt=dt;
    end
    up=U(f,i-1)+dt*V(f,i-1)+dt^2*(.5-beta)*A(f,i-1);
    vp=V(f,i-1)+dt*(1-gamma)*A(f,i-1);
    force=getforce(forceFcn,time(i),nd)-system.K(:,c)*bc.values;
    rhs=force(f)-C*vp-K*up; ai=zeros(numel(f),1);
    ai(order)=L'\(L\rhs(order));
    U(f,i)=up+beta*dt^2*ai; U(c,i)=bc.values;
    V(f,i)=vp+gamma*dt*ai; A(f,i)=ai;
    residual(i)=norm(M*ai+C*V(f,i)+K*U(f,i)-force(f),inf);
    energy(i)=.5*V(f,i)'*M*V(f,i)+.5*U(:,i)'*system.K*U(:,i);
end
result=struct('time',time,'U',U,'V',V,'A',A,'energy',energy, ...
    'freeResidual',residual,'rayleigh',rayleigh,'beta',beta,'gamma',gamma);
end

function F=getforce(fun,t,nd)
F=fun(t);
if isequal(size(F),[nd/3,3]), F=reshape(F',[],1); end
F=F(:); assert(numel(F)==nd && all(isfinite(F)),'fem3d:Load','Invalid force callback.');
end
