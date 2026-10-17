function [state,result] = initialize(nodes,velocity,flow,userOptions)
%INITIALIZE Bound rings and explicit material wake circulation state.
if nargin<4,userOptions=struct();end
o=uvlm.options(userOptions); g=uvlm.lattice(nodes); A=uvlm.influence(g,o.coreRadius);
TE=g.trailingEdge; direction=uvlm.air_velocity(flow,0,mean(TE,1)); direction=direction/norm(direction);
wake=struct('nodes',cat(1,reshape(TE,1,g.ny+1,3),reshape(TE+o.wakeLength*direction,1,g.ny+1,3)), ...
    'gamma',zeros(1,g.ny),'shedCirculation',zeros(0,g.ny));
gamma=zeros(g.nx*g.ny,1);
if o.steadyStart
    [V,R]=uvlm.wake_lattice(wake.nodes); coupled=A;
    for j=1:g.ny
        w=uvlm.ring_velocity(g.panels.control,V,R(j,:),1,o.coreRadius);
        coupled(:,j*g.nx)=coupled(:,j*g.nx)+sum(w.*g.panels.normal,2);
    end
    rhs=-sum((uvlm.air_velocity(flow,0,g.panels.control)-g.panels.Hcontrol*reshape(velocity,[],3)).*g.panels.normal,2);
    gamma=coupled\rhs; temp=reshape(gamma,g.nx,g.ny); wake.gamma=temp(end,:);
end
state=struct('time',0,'gamma',gamma,'wake',wake,'geometry',g,'influence',A,'options',o);
result=uvlm.loads(g,gamma,gamma,Inf,wake,flow,0,velocity,o);
result.boundaryResidual=0; result.kelvinDefect=0;
end
