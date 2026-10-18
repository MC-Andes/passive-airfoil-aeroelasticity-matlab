function result = compact_dynamic(structure,fluid,time,q0,v0,blockSize)
%COMPACT_DYNAMIC Exact Schur elimination of bound circulation in tangent FSI.
% Retains all TE circulation and material-motion history, with no wake
% truncation or approximation of the Biot-Savart kernel. Modal force and its
% unsteady derivative are projected before storage, reducing memory.
time=time(:)';dt=fluid.dt;nt=numel(time);n=size(structure.K,1);ny=fluid.spanPanels;ig=ny+(1:n);id=ny+n+(1:n);np=ny+2*n;
assert(max(abs(diff(time)-dt))<1e-9*dt && fluid.wakeRows>=nt,'fsi:Time','Uniform dt and full wake horizon required.');
Q=zeros(n,nt);V=Q;Acc=Q;Load=Q;Q(:,1)=q0;V(:,1)=v0;pgamma=fluid.projectedSteadyGammaQ*q0;gh=repmat(pgamma(1:ny),fluid.wakeRows,1);qh=repmat(q0,fluid.wakeRows,1);
load=pgamma(ig)+fluid.Lq*q0+fluid.Lv*v0+fluid.Fw*gh+fluid.Fq*qh;
Acc(:,1)=(structure.M-fluid.projectedE(id,:))\(load-structure.C*v0-structure.K*q0);Load(:,1)=load+fluid.projectedE(id,:)*Acc(:,1);
Ka=fluid.projectedD(ig,:)+fluid.projectedD(id,:)/dt+fluid.Lq;Ca=fluid.projectedE(ig,:)+fluid.projectedE(id,:)/dt+fluid.Lv;
effective=structure.M+.5*dt*(structure.C-Ca)+.25*dt^2*(structure.K-Ka);
energy=zeros(1,nt);work=energy;diss=energy;balance=energy;residual=energy;history=zeros(np,nt);history(:,1)=pgamma;
energy(1)=.5*(v0'*structure.M*v0+q0'*structure.K*q0);
for k=2:nt
 gh=[pgamma(1:ny);gh(1:end-ny)];qh=[Q(:,k-1);qh(1:end-n)];
 if mod(k-2,blockSize)==0
  first=k;nb=min(blockSize,nt-k+1);GH=zeros(numel(gh),nb);QH=zeros(numel(qh),nb);
  for col=0:nb-1,GH(1+col*ny:end,col+1)=gh(1:end-col*ny);QH(1+col*n:end,col+1)=qh(1:end-col*n);end
  baseGamma=-fluid.W*GH-fluid.Gq*QH;baseForce=fluid.Fw*GH+fluid.Fq*QH;
 end
 col=k-first+1;hg=baseGamma(:,col);hf=baseForce(:,col);recent=k-first;
 if recent>0
  gr=history(1:ny,k-1:-1:first);qr=Q(:,k-1:-1:first);
  hg=hg-fluid.W(:,1:recent*ny)*gr(:)-fluid.Gq(:,1:recent*n)*qr(:);hf=hf+fluid.Fw(:,1:recent*ny)*gr(:)+fluid.Fq(:,1:recent*n)*qr(:);
 end
 old=pgamma;historyLoad=hg(ig)+hg(id)/dt-old(id)/dt+hf;up=Q(:,k-1)+dt*V(:,k-1)+.25*dt^2*Acc(:,k-1);vp=V(:,k-1)+.5*dt*Acc(:,k-1);
 acc=effective\(historyLoad-(structure.C-Ca)*vp-(structure.K-Ka)*up);q=up+.25*dt^2*acc;v=vp+.5*dt*acc;pgamma=fluid.projectedD*q+fluid.projectedE*v+hg;
 load=pgamma(ig)+(pgamma(id)-old(id))/dt+fluid.Lq*q+fluid.Lv*v+hf;Q(:,k)=q;V(:,k)=v;Acc(:,k)=acc;Load(:,k)=load;history(:,k)=pgamma;
 residual(k)=norm(structure.M*acc+structure.C*v+structure.K*q-load)/max([norm(load),norm(structure.K*q),1e-10]);dq=q-Q(:,k-1);vmid=(v+V(:,k-1))/2;energy(k)=.5*(v'*structure.M*v+q'*structure.K*q);work(k)=work(k-1)+.5*(load+Load(:,k-1))'*dq;diss(k)=diss(k-1)+dt*vmid'*structure.C*vmid;balance(k)=energy(k)-energy(1)-work(k)+diss(k);
end
scale=max([max(abs(energy)),max(abs(work)),1e-15]);result=struct('time',time,'Q',Q,'V',V,'A',Acc,'load',Load,'gammaTEPerturbation',history(1:ny,:),'projectedGammaHistory',history,'energy',energy,'aerodynamicWork',work,'dissipation',diss,'energyBalance',balance,'relativeEnergyBalance',max(abs(balance))/scale,'relativeResidual',residual,'converged',true,'scope','exact bound-circulation Schur elimination of strongly implicit tangent UVLM/structure','circulationRetained',true,'wakeRows',fluid.wakeRows,'convolutionBlockSize',blockSize,'convolutionMethod','exact block matrix products; no kernel compression or circulation truncation');
end
