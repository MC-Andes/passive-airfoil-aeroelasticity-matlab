function result = transient(model,time,inputFcn,q0,v0)
%TRANSIENT Average-acceleration Newmark, prescribed time-dependent inputs.
% q0,v0 are modal coordinates. No pressure feedback or wake state.
time=time(:)'; n=size(model.K,1);
if nargin<4, q0=zeros(n,1); end
if nargin<5, v0=zeros(n,1); end
assert(numel(time)>1 && all(isfinite(time)) && all(diff(time)>0) && ...
    numel(q0)==n && numel(v0)==n && all(isfinite([q0(:);v0(:)])), ...
    'ael:Time','Invalid times or initial state.');
Q=zeros(n,numel(time)); V=Q; A=Q; Q(:,1)=q0(:); V(:,1)=v0(:);
M=model.M; K=model.K; C=model.C;
A(:,1)=M\(force(model,inputFcn,time(1))-C*V(:,1)-K*Q(:,1));
energy=zeros(1,numel(time)); residual=energy;
energy(1)=.5*(V(:,1)'*M*V(:,1)+Q(:,1)'*K*Q(:,1));
for k=2:numel(time)
    dt=time(k)-time(k-1); up=Q(:,k-1)+dt*V(:,k-1)+dt^2*.25*A(:,k-1);
    vp=V(:,k-1)+dt*.5*A(:,k-1); F=force(model,inputFcn,time(k));
    A(:,k)=(M+.5*dt*C+.25*dt^2*K)\(F-C*vp-K*up);
    Q(:,k)=up+.25*dt^2*A(:,k); V(:,k)=vp+.5*dt*A(:,k);
    energy(k)=.5*(V(:,k)'*M*V(:,k)+Q(:,k)'*K*Q(:,k));
    residual(k)=norm(M*A(:,k)+C*V(:,k)+K*Q(:,k)-F,inf);
end
result=struct('time',time,'Q',Q,'V',V,'A',A,'Y',model.O*Q, ...
    'energy',energy,'residual',residual,'scope','prescribed-input structural response');
end

function F=force(model,fun,t)
p=fun(t); p=p(:);
assert(numel(p)==size(model.B,2) && all(isfinite(p)), ...
    'ael:Load','Input callback must return one amplitude per input channel.');
F=model.B*p;
end
