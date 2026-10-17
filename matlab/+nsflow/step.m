function [next,info] = step(g,state,dt,physical)
%STEP Projected explicit midpoint RK2, with optional Brinkman wall penalization.
% No-slip is an APPROXIMATION when an immersed mask is used; eta and geometry
% resolution must be refined. Without a mask this is a second-order MAC solver.
nu=physical.nu;rho=physical.rho;d=g.dimension;U=state.U;
assert(nu>0 && rho>0 && dt>0,'nsflow:Physical','Positive viscosity, density and dt required.');
courant=0;for i=1:d,courant=courant+dt*max(abs(U{i}),[],'all')/g.h(i);end
assert(courant<=.8 && dt*nu*sum(1./g.h.^2)<=.45,'nsflow:Time','Refine dt for advection/diffusion stability.');
force=external(state.time,U);R=nsflow.rhs(g,U,nu,force);predict=cell(1,d);
for i=1:d,predict{i}=U{i}+.5*dt*R{i};end
coupled=isfield(physical,'mask') && (~isfield(physical,'penaltyProjection') || strcmp(physical.penaltyProjection,'coupled'));
if coupled
 [predict,~,~]=nsflow.project_penalized(g,predict,state.time+.5*dt,.5*dt,physical);
else
 [predict,~,~]=penalize(predict,state.time+.5*dt,.5*dt);[predict,~,~]=nsflow.project(g,predict,.5*dt,rho);
end
midForce=external(state.time+.5*dt,predict);R=nsflow.rhs(g,predict,nu,midForce);candidate=cell(1,d);
for i=1:d,candidate{i}=U{i}+dt*R{i};end
if coupled
 [candidate,p,projection]=nsflow.project_penalized(g,candidate,state.time+dt,dt,physical);bodyForce=projection.bodyForce;wallWork=projection.wallWork;divergence=projection.maximumDivergence;
else
 [candidate,bodyForce,wallWork]=penalize(candidate,state.time+dt,dt);[candidate,p,divergence]=nsflow.project(g,candidate,dt,rho);
end
next=struct('U',{candidate},'p',p,'time',state.time+dt);
kinetic=0;for i=1:d,kinetic=kinetic+.5*rho*g.cellVolume*sum(candidate{i}.^2,'all');end
reservoir=zeros(1,d);momentumChange=zeros(1,d);
for i=1:d,reservoir(i)=rho*g.cellVolume*sum(midForce{i},'all');momentumChange(i)=rho*g.cellVolume*sum(candidate{i}-U{i},'all');end
momentumDefect=momentumChange-dt*(reservoir-bodyForce);
info=struct('momentumBalance',momentumDefect,'reservoirForce',reservoir,'maximumDivergence',divergence,'kineticEnergy',kinetic,'bodyForce',bodyForce, ...
    'wallWork',wallWork,'courant',courant,'pressureGauge','zero periodic mean');
    function F=external(t,vel)
        F=cell(1,d);for a=1:d,F{a}=zeros(g.count);end
        if isfield(physical,'acceleration'),for a=1:d,F{a}=F{a}+physical.acceleration(a);end,end
        if isfield(physical,'forcing'),F=physical.forcing(g,t,vel);end
        if isfield(physical,'fringe')
            for a=1:d,rate=physical.fringe(g.velocityCoordinates(a,:));F{a}=F{a}+rate.*(physical.freestream(a)-vel{a});end
        end
    end
    function [vel,bodyForce,work]=penalize(vel,t,h)
        bodyForce=zeros(1,d);work=0;if ~isfield(physical,'mask'),return;end
        assert(isfield(physical,'penaltyTime') && physical.penaltyTime>0,'nsflow:Wall','Declare the Brinkman penalty time.');
        for a=1:d
            mask=physical.mask(g.velocityCoordinates(a,:),t);wall=zeros(g.count);
            if isfield(physical,'wallVelocity'),wall=physical.wallVelocity(g.velocityCoordinates(a,:),t,a);end
            before=vel{a};vel{a}=(before+h*mask.*wall/physical.penaltyTime)./(1+h*mask/physical.penaltyTime);
            densityForce=rho*(vel{a}-before)/h;bodyForce(a)=-g.cellVolume*sum(densityForce,'all');
            work=work+g.cellVolume*sum(densityForce.*wall,'all');
        end
    end
end
