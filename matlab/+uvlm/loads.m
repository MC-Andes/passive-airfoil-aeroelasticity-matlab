function r = loads(g,gamma,previousGamma,dt,wake,flow,t,velocity,options)
%LOADS Bound-vorticity Joukowski forces plus unsteady Bernoulli contribution.
% Net bound edges exclude the shed TE closing segment. No skin Cp claim.
% Force points have explicit material interpolation and adjoint nodal loads.
G=reshape(gamma,g.nx,g.ny); spanStrength=G-[zeros(1,g.ny);G(1:end-1,:)];
chordStrength=[zeros(g.nx,1),G]-[G,zeros(g.nx,1)];
[wV,wR]=uvlm.wake_lattice(wake.nodes); U=reshape(velocity,[],3);
Pspan=g.Hspan*reshape(g.nodes,[],3); Pchord=g.Hchord*reshape(g.nodes,[],3);
allPoints=[Pspan;Pchord]; ambient=uvlm.air_velocity(flow,t,allPoints);
induced=uvlm.ring_velocity(allPoints,g.ringVertices,g.ring,gamma,options.coreRadius) ...
    +uvlm.ring_velocity(allPoints,wV,wR,wake.gamma(:),options.coreRadius);
ns=size(Pspan,1); relative=ambient+induced-[g.Hspan*U;g.Hchord*U];
Fspan=flow.rho*spanStrength(:).*cross(relative(1:ns,:),g.spanB-g.spanA,2);
Fchord=flow.rho*chordStrength(:).*cross(relative(ns+1:end,:),g.chordB-g.chordA,2);
dG=zeros(size(gamma)); if isfinite(dt), dG=(gamma-previousGamma)/dt; end
Funsteady=flow.rho*dG.*g.panels.area.*g.panels.normal;
nodal=g.Hspan'*Fspan+g.Hchord'*Fchord+g.Hcent'*Funsteady;
Fpanel=Fspan+Funsteady; chord3=reshape(Fchord,g.nx,g.ny+1,3);
for j=1:g.ny+1
    if j==1, ids=1:g.nx; coefficient=1;
    elseif j==g.ny+1, ids=(g.ny-1)*g.nx+(1:g.nx); coefficient=1;
    else, ids=(j-2)*g.nx+(1:2*g.nx); coefficient=.5; end
    value=reshape(chord3(:,j,:),g.nx,3);
    if numel(ids)==2*g.nx, value=[value;value]; end
    Fpanel(ids,:)=Fpanel(ids,:)+coefficient*value;
end
X=reshape(g.nodes,[],3); total=sum(nodal,1); ref=[0,0,0]; if isfield(flow,'referencePoint'),ref=flow.referencePoint;end
moment=sum(cross(X-ref,nodal,2),1); direction=uvlm.air_velocity(rmfield_if(flow,'gust'),t,[0,0,0]); direction=direction/norm(direction);
liftDirection=cross(direction,[0,1,0]); liftDirection=liftDirection/norm(liftDirection);
area=sum(g.panels.area); if isfield(flow,'referenceArea'), area=flow.referenceArea; end
chord=max(X(:,1))-min(X(:,1)); if isfield(flow,'referenceChord'),chord=flow.referenceChord;end
qS=.5*flow.rho*flow.speed^2*area;
r=struct('gamma',gamma,'gammaDot',dG,'nodalForces',nodal,'panelForces',Fpanel, ...
    'totalForce',total,'totalMoment',moment,'lift',dot(total,liftDirection),'CL',dot(total,liftDirection)/qS, ...
    'CDi',dot(total,direction)/qS,'Cm',moment(2)/(qS*chord), ...
    'unsteadyForce',sum(Funsteady,1),'pressureJump',sum(Fpanel.*g.panels.normal,2)./g.panels.area, ...
    'fidelity','inviscid UVLM; Joukowski bound edges plus unsteady Bernoulli; not viscous CFD');
end
function s=rmfield_if(s,name)
if isfield(s,name),s=rmfield(s,name);end
end
