function model = linear_model(nodes,Phi,flow,dt,steps,settings)
%LINEAR_MODEL Tangent UVLM about a rigid/deformed equilibrium with frozen wake.
% Gamma history AND material trailing-edge displacement history are retained.
% Every actually shed row fits in the buffer (steps+1); no emitted circulation
% is truncated. Constant prehistory beyond the horizon is a finite-wake limit.
if nargin<6,settings=struct();end
wakeDerivative='analytic';if isfield(settings,'wakeDerivative'),wakeDerivative=settings.wakeDerivative;settings=rmfield(settings,'wakeDerivative');end
compact=false;if isfield(settings,'compactTangent'),compact=settings.compactTangent;settings=rmfield(settings,'compactTangent');end
batchSize=64;if isfield(settings,'ageBatchSize'),batchSize=settings.ageBatchSize;settings=rmfield(settings,'ageBatchSize');end
o=uvlm.options(settings);assert(strcmp(o.convection,'frozen'),'uvlm:Linear','Tangent model requires frozen ambient convection.');
[state,base]=uvlm.initialize(nodes,zeros(numel(nodes)/3,3),flow,o);g=state.geometry;
np=g.nx*g.ny;nm=size(Phi,2);nn=numel(nodes)/3;N=steps+1;X=reshape(nodes,[],3);
A=state.influence;Bq=zeros(np,nm);Bv=Bq;Lq=zeros(nm);Lv=Lq;Lg=zeros(nm,np);Ld=Lg;
TEindices=reshape(1:nn,g.nx+1,g.ny+1);TEindices=TEindices(end,:);
TEdofs=reshape([3*TEindices-2;3*TEindices-1;3*TEindices],[],1);TEPhi=Phi(TEdofs,:);
TEextract=sparse(1:g.ny,g.nx:g.nx:np,1,g.ny,np);gamma0=state.gamma;G0=reshape(gamma0,g.nx,g.ny);
span0=G0-[zeros(1,g.ny);G0(1:end-1,:)];chord0=[zeros(g.nx,1),G0]-[G0,zeros(g.nx,1)];
Pspan=g.Hspan*X;Pchord=g.Hchord*X;points=[g.panels.control;Pspan;Pchord];nc=np;ns=size(Pspan,1);
XiSpan=kron(g.Hspan,speye(3))*Phi;XiChord=kron(g.Hchord,speye(3))*Phi;
geometryStep=zeros(nm,1);steadyGammaQ=zeros(np,nm);
for mode=1:nm
    h=1e-5*flow.referenceChord/max(max(abs(Phi(:,mode))),1e-20);geometryStep(mode)=h;
    xp=reshape(X+reshape(h*Phi(:,mode),3,[])',size(nodes));xm=reshape(X-reshape(h*Phi(:,mode),3,[])',size(nodes));
    [washp,fp]=shape(xp);[washm,fm]=shape(xm);Bq(:,mode)=(washp-washm)/(2*h);Lq(:,mode)=Phi'*(fp-fm)/(2*h);
    vp=reshape(Phi(:,mode),3,[])';Bv(:,mode)=sum((g.panels.Hcontrol*vp).*g.panels.normal,2);
    ap=uvlm.loads(g,gamma0,gamma0,Inf,state.wake,flow,0,vp,o);am=uvlm.loads(g,gamma0,gamma0,Inf,state.wake,flow,0,-vp,o);
    Lv(:,mode)=Phi'*(reshape(ap.nodalForces',[],1)-reshape(am.nodalForces',[],1))/2;
    [sp,~]=uvlm.initialize(xp,zeros(nn,3),flow,o);[sm,~]=uvlm.initialize(xm,zeros(nn,3),flow,o);
    steadyGammaQ(:,mode)=(sp.gamma-sm.gamma)/(2*h);
end
% Exact derivative of the quadratic Joukowski load: direct bound-edge
% strength plus baseline circulation times induced-velocity derivative.
% This avoids repeated complete load evaluations for every circulation DOF.
forcePoints=[Pspan;Pchord];[wv,wr]=uvlm.wake_lattice(state.wake.nodes);
relative=uvlm.air_velocity(flow,0,forcePoints)+uvlm.ring_velocity(forcePoints,g.ringVertices,g.ring,gamma0,o.coreRadius)+uvlm.ring_velocity(forcePoints,wv,wr,state.wake.gamma(:),o.coreRadius);
fs=flow.rho*cross(relative(1:ns,:),g.spanB-g.spanA,2);fc=flow.rho*cross(relative(ns+1:end,:),g.chordB-g.chordA,2);
Qs=reshape(sum(permute(reshape(XiSpan,3,ns,nm),[2,1,3]).*fs,2),ns,nm)';nch=size(Pchord,1);Qc=reshape(sum(permute(reshape(XiChord,3,nch,nm),[2,1,3]).*fc,2),nch,nm)';
ids=(1:np)';next=ids+1;take=mod(ids,g.nx)~=0;Sspan=speye(np)+sparse(next(take),ids(take),-1,np,np);Schord=sparse([ids;ids+g.nx],[ids;ids],[-ones(np,1);ones(np,1)],nch,np);
Lg=Qs*Sspan+Qc*Schord;
for first=1:64:np
    panels=first:min(first+63,np);rings=g.ring(panels,:);induced=zeros(size(forcePoints,1),3,numel(panels));
    for edge=1:4,induced=induced+uvlm.segment_influence(forcePoints,g.ringVertices(rings(:,edge),:),g.ringVertices(rings(:,mod(edge,4)+1),:),o.coreRadius);end
    Lg(:,panels)=Lg(:,panels)+project_batch(induced);
end
XiCent=kron(g.Hcent,speye(3))*Phi;direct=flow.rho*g.panels.area.*g.panels.normal;Ld=reshape(sum(permute(reshape(XiCent,3,np,nm),[2,1,3]).*direct,2),np,nm)';
projection=eye(np);if compact,projection=[full(TEextract);Lg;Ld]/A;end
nh=size(projection,1);W=zeros(nh,g.ny,N);Fw=zeros(nm,g.ny,N);Gq=zeros(nh,nm,N);Fq=zeros(nm,nm,N);
batchSize=max(1,min(batchSize,floor(120000/size(points,1))));
ambient=uvlm.air_velocity(flow,0,[0,0,0]);TE=g.trailingEdge;lambda=[0,G0(end,:)]-[G0(end,:),0];
for first=1:batchSize:N
    if first==1 || mod(first-1,200*batchSize)==0,fprintf('UVLM wake kernel ages %d/%d (%s)\n',first,N,wakeDerivative);end
    ages=first:min(first+batchSize-1,N);nb=numel(ages);ny1=g.ny+1;
    for strip=1:g.ny
        a=TE(strip,:)+(ages'-1)*dt*ambient;b=TE(strip+1,:)+(ages'-1)*dt*ambient;c=b+dt*ambient;d=a+dt*ambient;
        induced=uvlm.segment_influence(points,a,b,o.coreRadius)+uvlm.segment_influence(points,b,c,o.coreRadius)+uvlm.segment_influence(points,c,d,o.coreRadius)+uvlm.segment_influence(points,d,a,o.coreRadius);
        wash=reshape(sum(induced(1:nc,:,:).*g.panels.normal,2),np,nb);W(:,strip,ages)=reshape(projection*wash,nh,1,nb);
        Fw(:,strip,ages)=reshape(project_batch(induced(nc+1:end,:,:)),nm,1,nb);
    end
    if strcmp(wakeDerivative,'analytic')
        for strip=1:ny1
            if lambda(strip)==0,continue;end
            center=TE(strip,:)+ages'*dt*ambient;derivative=uvlm.wake_vertex_gradient(points,center,dt*ambient,o.coreRadius);
            for axis=1:3
                induced=derivative{axis};wash=reshape(sum(induced(1:nc,:,:).*g.panels.normal,2),np,nb);dw=projection*wash;df=project_batch(induced(nc+1:end,:,:));weights=TEPhi(3*(strip-1)+axis,:);
                Gq(:,:,ages)=Gq(:,:,ages)+lambda(strip)*reshape(dw,nh,1,nb).*reshape(weights,1,nm,1);Fq(:,:,ages)=Fq(:,:,ages)+lambda(strip)*reshape(df,nm,1,nb).*reshape(weights,1,nm,1);
            end
        end
    elseif strcmp(wakeDerivative,'complex')
        % The regularized kernel is analytic in endpoint coordinates away
        % from its declared cutoff. Complex-step differentiation avoids
        % subtraction and reuses three coordinate derivatives for all modes.
        h=1e-20*flow.referenceChord;
        for strip=1:ny1
            if lambda(strip)==0,continue;end
            center=TE(strip,:)+ages'*dt*ambient;previous=center-dt*ambient;following=center+dt*ambient;
            for axis=1:3
                cp=complex(center);cp(:,axis)=cp(:,axis)+1i*h;
                induced=imag(uvlm.segment_influence(points,previous,cp,o.coreRadius)+uvlm.segment_influence(points,cp,following,o.coreRadius))/h;
                wash=reshape(sum(induced(1:nc,:,:).*g.panels.normal,2),np,nb);dw=projection*wash;df=project_batch(induced(nc+1:end,:,:));weights=TEPhi(3*(strip-1)+axis,:);
                Gq(:,:,ages)=Gq(:,:,ages)+lambda(strip)*reshape(dw,nh,1,nb).*reshape(weights,1,nm,1);
                Fq(:,:,ages)=Fq(:,:,ages)+lambda(strip)*reshape(df,nm,1,nb).*reshape(weights,1,nm,1);
            end
        end
    else
    for mode=1:nm
        h=geometryStep(mode);move=reshape(TEPhi(:,mode),3,[])';vplus=zeros(size(points,1),3,nb);vminus=vplus;
        for strip=1:ny1
            if lambda(strip)==0,continue;end
            center=TE(strip,:)+ages'*dt*ambient;previous=center-dt*ambient;following=center+dt*ambient;cp=center+h*move(strip,:);cm=center-h*move(strip,:);
            vplus=vplus+lambda(strip)*(uvlm.segment_influence(points,previous,cp,o.coreRadius)+uvlm.segment_influence(points,cp,following,o.coreRadius));
            vminus=vminus+lambda(strip)*(uvlm.segment_influence(points,previous,cm,o.coreRadius)+uvlm.segment_influence(points,cm,following,o.coreRadius));
        end
        induced=(vplus-vminus)/(2*h);wash=reshape(sum(induced(1:nc,:,:).*g.panels.normal,2),np,nb);Gq(:,mode,ages)=reshape(projection*wash,nh,1,nb);Fq(:,mode,ages)=reshape(project_batch(induced(nc+1:end,:,:)),nm,1,nb);
    end
    end
end
model=struct('A',A,'D',A\(-Bq),'E',A\Bv,'W',reshape(W,nh,[]),'Fw',reshape(Fw,nm,[]), ...
    'Gq',reshape(Gq,nh,[]),'Fq',reshape(Fq,nm,[]),'Lg',Lg,'Ld',Ld,'Lq',Lq,'Lv',Lv, ...
    'TEextract',TEextract,'steadyGammaQ',steadyGammaQ,'gamma0',gamma0,'load0',Phi'*reshape(base.nodalForces',[],1), ...
    'base',base,'state0',state,'dt',dt,'wakeRows',N,'modes',nm,'spanPanels',g.ny, ...
    'wakeDerivative',wakeDerivative,'compactTangent',compact,'kernelBatchSize',batchSize,'PhiSurface',Phi,'X0',nodes,'flow',flow,'settings',o, ...
    'scope','tangent of own UVLM; frozen wake, all emitted circulation retained; small perturbations');
if compact
    P=[full(TEextract);Lg;Ld];model.projectedD=P*model.D;model.projectedE=P*model.E;model.projectedSteadyGammaQ=P*steadyGammaQ;
    model.projectionDefinition='[TE extraction; force circulation derivative Lg; added-mass derivative Ld] A^-1; exact Schur elimination';
end
    function [wash,force]=shape(x)
        gx=uvlm.lattice(x);wake=state.wake;
        near=reshape(g.trailingEdge+dt*uvlm.air_velocity(flow,0,[0,0,0]),1,g.ny+1,3);
        wake.nodes=cat(1,reshape(gx.trailingEdge,1,g.ny+1,3),near,state.wake.nodes(end,:,:));
        wake.gamma=[state.wake.gamma;state.wake.gamma];
        [wv,wr]=uvlm.wake_lattice(wake.nodes);
        induced=uvlm.ring_velocity(gx.panels.control,gx.ringVertices,gx.ring,gamma0,o.coreRadius)+ ...
            uvlm.ring_velocity(gx.panels.control,wv,wr,wake.gamma(:),o.coreRadius);
        wash=sum((uvlm.air_velocity(flow,0,gx.panels.control)+induced).*gx.panels.normal,2);
        a=uvlm.loads(gx,gamma0,gamma0,Inf,wake,flow,0,zeros(nn,3),o);force=reshape(a.nodalForces',[],1);
    end
    function Q=project_batch(induced)
        fs=flow.rho*span0(:).*cross(induced(1:ns,:,:),repmat(g.spanB-g.spanA,1,1,size(induced,3)),2);
        fc=flow.rho*chord0(:).*cross(induced(ns+1:end,:,:),repmat(g.chordB-g.chordA,1,1,size(induced,3)),2);
        Q=XiSpan'*reshape(permute(fs,[2,1,3]),3*ns,[])+XiChord'*reshape(permute(fc,[2,1,3]),3*size(Pchord,1),[]);
    end
end
