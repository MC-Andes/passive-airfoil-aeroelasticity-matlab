function model = refine_time_model(parent)
%REFINE_TIME_MODEL Exact half-step kernels using geometric vortex subdivision.
% Each coarse ring equals its two half rings. A coarse shared vertex moves
% its two intermediate half vertices by half its displacement. Therefore
% G_old(a)=G_new(2a)+(G_new(2a-1)+G_new(2a+1))/2. Compute odd kernels directly
% and recover even kernels by these exact identities, without approximating
% the new circulation or material-motion history. Geometry/flow stay fixed.
assert(parent.compactTangent && ~isfield(parent.flow,'gust'),'uvlm:Refine','Requires a compact steady-ambient frozen-wake model.');
model=parent;g=parent.state0.geometry;flow=parent.flow;dt=parent.dt/2;Nold=parent.wakeRows;N=2*(Nold-1);nm=parent.modes;ny=parent.spanPanels;np=g.nx*g.ny;nh=ny+2*nm;
X=reshape(parent.X0,[],3);Phi=parent.PhiSurface;Pspan=g.Hspan*X;Pchord=g.Hchord*X;points=[g.panels.control;Pspan;Pchord];ns=size(Pspan,1);nch=size(Pchord,1);
XiSpan=kron(g.Hspan,speye(3))*Phi;XiChord=kron(g.Hchord,speye(3))*Phi;G0=reshape(parent.gamma0,g.nx,g.ny);span0=G0-[zeros(1,ny);G0(1:end-1,:)];chord0=[zeros(g.nx,1),G0]-[G0,zeros(g.nx,1)];
projection=[full(parent.TEextract);parent.Lg;parent.Ld]/parent.A;
ids=reshape(1:size(X,1),g.nx+1,ny+1);te=ids(end,:);TEPhi=Phi(reshape([3*te-2;3*te-1;3*te],[],1),:);TE=g.trailingEdge;ambient=uvlm.air_velocity(flow,0,[0,0,0]);lambda=[0,G0(end,:)]-[G0(end,:),0];core=parent.settings.coreRadius;
Wo=zeros(nh,ny,Nold);Fo=zeros(nm,ny,Nold);Go=zeros(nh,nm,Nold);Qo=zeros(nm,nm,Nold);rawFirst=zeros(np,nm);batch=max(1,min(64,floor(120000/size(points,1))));
for first=1:batch:Nold
 if first==1 || mod(first-1,200*batch)==0,fprintf('UVLM exact subdivision odd kernels %d/%d\n',first,Nold);end
 slots=first:min(first+batch-1,Nold);ages=2*slots-1;nb=numel(slots);
 for strip=1:ny
  a=TE(strip,:)+(ages'-1)*dt*ambient;b=TE(strip+1,:)+(ages'-1)*dt*ambient;c=b+dt*ambient;d=a+dt*ambient;
  V=uvlm.segment_influence(points,a,b,core)+uvlm.segment_influence(points,b,c,core)+uvlm.segment_influence(points,c,d,core)+uvlm.segment_influence(points,d,a,core);
  wash=reshape(sum(V(1:np,:,:).*g.panels.normal,2),np,nb);Wo(:,strip,slots)=reshape(projection*wash,nh,1,nb);Fo(:,strip,slots)=reshape(force(V(np+1:end,:,:)),nm,1,nb);
 end
 for strip=1:ny+1
  if lambda(strip)==0,continue;end
  derivative=uvlm.wake_vertex_gradient(points,TE(strip,:)+ages'*dt*ambient,dt*ambient,core);
  for axis=1:3
   V=derivative{axis};wash=reshape(sum(V(1:np,:,:).*g.panels.normal,2),np,nb);dw=projection*wash;df=force(V(np+1:end,:,:));weights=TEPhi(3*(strip-1)+axis,:);
   Go(:,:,slots)=Go(:,:,slots)+lambda(strip)*reshape(dw,nh,1,nb).*reshape(weights,1,nm,1);Qo(:,:,slots)=Qo(:,:,slots)+lambda(strip)*reshape(df,nm,1,nb).*reshape(weights,1,nm,1);
   if first==1,rawFirst=rawFirst+lambda(strip)*wash(:,1)*weights;end
  end
 end
end
W=zeros(nh,ny,N);Fw=zeros(nm,ny,N);Gq=zeros(nh,nm,N);Fq=zeros(nm,nm,N);
oldW=reshape(parent.W,nh,ny,Nold);oldF=reshape(parent.Fw,nm,ny,Nold);oldG=reshape(parent.Gq,nh,nm,Nold);oldQ=reshape(parent.Fq,nm,nm,Nold);
W(:,:,1:2:end)=Wo(:,:,1:end-1);W(:,:,2:2:end)=oldW(:,:,1:end-1)-Wo(:,:,1:end-1);
Fw(:,:,1:2:end)=Fo(:,:,1:end-1);Fw(:,:,2:2:end)=oldF(:,:,1:end-1)-Fo(:,:,1:end-1);
Gq(:,:,1:2:end)=Go(:,:,1:end-1);Gq(:,:,2:2:end)=oldG(:,:,1:end-1)-(Go(:,:,1:end-1)+Go(:,:,2:end))/2;
Fq(:,:,1:2:end)=Qo(:,:,1:end-1);Fq(:,:,2:2:end)=oldQ(:,:,1:end-1)-(Qo(:,:,1:end-1)+Qo(:,:,2:end))/2;
model.W=reshape(W,nh,[]);model.Fw=reshape(Fw,nm,[]);model.Gq=reshape(Gq,nh,[]);model.Fq=reshape(Fq,nm,[]);
model.D=parent.D+.5*(parent.A\rawFirst);model.projectedD=parent.projectedD+.5*Go(:,:,1);model.Lq=parent.Lq-.5*Qo(:,:,1);model.dt=dt;model.wakeRows=N;model.kernelBatchSize=batch;
model.refinement=struct('parentDt',parent.dt,'parentWakeRows',Nold,'method','exact spatial subdivision identities; odd kernels evaluated directly; independent fine histories fully retained');
 function Q=force(V)
  fs=flow.rho*span0(:).*cross(V(1:ns,:,:),repmat(g.spanB-g.spanA,1,1,size(V,3)),2);fc=flow.rho*chord0(:).*cross(V(ns+1:end,:,:),repmat(g.chordB-g.chordA,1,1,size(V,3)),2);
  Q=XiSpan'*reshape(permute(fs,[2,1,3]),3*ns,[])+XiChord'*reshape(permute(fc,[2,1,3]),3*nch,[]);
 end
end
