function [next,result] = step(state,nodes,velocity,flow,dt)
%STEP Material wake advection, Kutta shedding and current bound circulation.
% Wake is copied from the OLD state: safe for FSI subiterations. No truncation.
assert(isscalar(dt) && dt>0 && isfinite(dt),'uvlm:Time','Positive finite dt required.');
o=state.options;
chords=vecnorm(diff(nodes,1,1),2,3);
assert(dt*(flow.speed+max(vecnorm(reshape(velocity,[],3),2,2)))/min(chords,[],'all')<=o.courantLimit, ...
    'uvlm:Courant','Refine dt: wake travel exceeds the configured minimum-panel chord ratio.');
t=state.time; old=state.geometry; wake=state.wake;
assert(size(wake.gamma,1)<o.maximumWakeRows,'uvlm:WakeBudget', ...
    'Wake budget exhausted; circulation will not be silently truncated.');
P=reshape(wake.nodes,[],3); speed=convection(P,t);
if o.convectionOrder==2, speed=(speed+convection(P+dt*speed,t+dt))/2; end
convected=reshape(P+dt*speed,size(wake.nodes));
g=uvlm.lattice(nodes); trailing=reshape(state.gamma,old.nx,old.ny); shed=trailing(end,:);
wake.nodes=cat(1,reshape(g.trailingEdge,1,g.ny+1,3),convected); wake.gamma=[shed;wake.gamma];
wake.shedCirculation=[shed;wake.shedCirculation];
[V,R]=uvlm.wake_lattice(wake.nodes);
induced=uvlm.ring_velocity(g.panels.control,V,R,wake.gamma(:),o.coreRadius);
if isequal(nodes,old.nodes), A=state.influence; else,A=uvlm.influence(g,o.coreRadius);end
rhs=-sum((uvlm.air_velocity(flow,t+dt,g.panels.control)+induced-g.panels.Hcontrol*reshape(velocity,[],3)).*g.panels.normal,2);
gamma=A\rhs;
next=struct('time',t+dt,'gamma',gamma,'wake',wake,'geometry',g,'influence',A,'options',o);
result=uvlm.loads(g,gamma,state.gamma,dt,wake,flow,t+dt,velocity,o);
result.boundaryResidual=norm(A*gamma-rhs,inf)/max(flow.speed,1);
result.kelvinDefect=max(abs(wake.gamma(1,:)-shed));
assert(isequal(wake.gamma(2:end,:),state.wake.gamma),'uvlm:Kelvin','Shed circulation was modified.');
    function speed=convection(points,time)
        speed=uvlm.air_velocity(flow,time,points);
        if strcmp(o.convection,'free')
            [v0,r0]=uvlm.wake_lattice(state.wake.nodes);
            speed=speed+uvlm.ring_velocity(points,old.ringVertices,old.ring,state.gamma,o.coreRadius) ...
                +uvlm.ring_velocity(points,v0,r0,state.wake.gamma(:),o.coreRadius);
        end
    end
end
