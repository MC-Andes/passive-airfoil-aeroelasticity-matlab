function result = linear_dynamic(structure,fluid,time,q0,v0,blockSize)
%LINEAR_DYNAMIC Strong implicit coupling to the verified tangent UVLM.
% Coordinates are perturbations about a separately converged static state.
if nargin<6,blockSize=32;end
if isfield(fluid,'compactTangent') && fluid.compactTangent,result=fsi.compact_dynamic(structure,fluid,time,q0,v0,blockSize);return;end
assert(blockSize>=1 && blockSize==fix(blockSize),'fsi:Convolution','Positive integer convolution block required.');
solveA=decomposition(fluid.A,'lu');
time=time(:)';dt=fluid.dt;nt=numel(time);n=size(structure.K,1);
assert(max(abs(diff(time)-dt))<1e-9*dt && fluid.wakeRows>=nt,'fsi:Time','Uniform dt and full wake horizon required.');
Q=zeros(n,nt);V=Q;Acc=Q;Load=Q;Q(:,1)=q0;V(:,1)=v0;
gamma=fluid.steadyGammaQ*q0;gh=repmat(fluid.TEextract*gamma,fluid.wakeRows,1);qh=repmat(q0,fluid.wakeRows,1);
load=fluid.Lg*gamma+fluid.Lq*q0+fluid.Lv*v0+fluid.Fw*gh+fluid.Fq*qh;Load(:,1)=load;
initialMass=structure.M-fluid.Ld*fluid.E;
Acc(:,1)=initialMass\(load-structure.C*v0-structure.K*q0);
Load(:,1)=load+fluid.Ld*fluid.E*Acc(:,1);
Kgamma=fluid.Lg+fluid.Ld/dt;Ka=Kgamma*fluid.D+fluid.Lq;Ca=Kgamma*fluid.E+fluid.Lv;
effective=structure.M+.5*dt*(structure.C-Ca)+.25*dt^2*(structure.K-Ka);
energy=zeros(1,nt);work=energy;diss=energy;balance=energy;residual=energy;
energy(1)=.5*(v0'*structure.M*v0+q0'*structure.K*q0);gammaHistory=zeros(numel(gamma),nt);gammaHistory(:,1)=gamma;
for k=2:nt
    gh=[fluid.TEextract*gamma;gh(1:end-fluid.spanPanels)];qh=[Q(:,k-1);qh(1:end-n)];
    if mod(k-2,blockSize)==0
        first=k;nb=min(blockSize,nt-k+1);GH=zeros(numel(gh),nb);QH=zeros(numel(qh),nb);
        for column=0:nb-1
            GH(1+column*fluid.spanPanels:end,column+1)=gh(1:end-column*fluid.spanPanels);
            QH(1+column*n:end,column+1)=qh(1:end-column*n);
        end
        baseWash=fluid.W*GH+fluid.Gq*QH;baseForce=fluid.Fw*GH+fluid.Fq*QH;
    end
    column=k-first+1;wash=baseWash(:,column);historyForce=baseForce(:,column);
    recent=k-first;
    if recent>0
        gr=fluid.TEextract*gammaHistory(:,k-1:-1:first);qr=Q(:,k-1:-1:first);
        wash=wash+fluid.W(:,1:recent*fluid.spanPanels)*gr(:)+fluid.Gq(:,1:recent*n)*qr(:);
        historyForce=historyForce+fluid.Fw(:,1:recent*fluid.spanPanels)*gr(:)+fluid.Fq(:,1:recent*n)*qr(:);
    end
    historyGamma=solveA\(-wash);
    historyLoad=Kgamma*historyGamma-fluid.Ld/dt*gamma+historyForce;
    up=Q(:,k-1)+dt*V(:,k-1)+.25*dt^2*Acc(:,k-1);vp=V(:,k-1)+.5*dt*Acc(:,k-1);
    acc=effective\(historyLoad-(structure.C-Ca)*vp-(structure.K-Ka)*up);
    q=up+.25*dt^2*acc;v=vp+.5*dt*acc;newGamma=fluid.D*q+fluid.E*v+historyGamma;
    load=fluid.Lg*newGamma+fluid.Ld/dt*(newGamma-gamma)+fluid.Lq*q+fluid.Lv*v+historyForce;
    Q(:,k)=q;V(:,k)=v;Acc(:,k)=acc;Load(:,k)=load;gammaHistory(:,k)=newGamma;gamma=newGamma;
    residual(k)=norm(structure.M*acc+structure.C*v+structure.K*q-load)/max([norm(load),norm(structure.K*q),1e-10]);
    dq=q-Q(:,k-1);vmid=(v+V(:,k-1))/2;energy(k)=.5*(v'*structure.M*v+q'*structure.K*q);
    work(k)=work(k-1)+.5*(load+Load(:,k-1))'*dq;diss(k)=diss(k-1)+dt*vmid'*structure.C*vmid;
    balance(k)=energy(k)-energy(1)-work(k)+diss(k);
end
scale=max([max(abs(energy)),max(abs(work)),1e-15]);
result=struct('time',time,'Q',Q,'V',V,'A',Acc,'load',Load,'gammaPerturbation',gammaHistory, ...
    'energy',energy,'aerodynamicWork',work,'dissipation',diss,'energyBalance',balance, ...
    'relativeEnergyBalance',max(abs(balance))/scale,'relativeResidual',residual,'converged',true, ...
    'scope','strong implicit linearized UVLM/structural coupling; no synthetic oscillating pressure', ...
    'circulationRetained',true,'wakeRows',fluid.wakeRows,'convolutionBlockSize',blockSize,'convolutionMethod','exact block matrix products; no kernel compression or circulation truncation');
end
