function response = simulate(model,flowInput,time,q0,v0,overrides)
%SIMULATE Implicit Newmark integration with two-way quasi-steady VLM loads.
% Constant reduced M,C; nonlinear restoring force and updated surface.
% Use only when aerodynamic memory can be neglected; no UVLM/added mass.
if nargin<4 || isempty(q0), q0=zeros(model.ndof,1); end
if nargin<5 || isempty(v0), v0=zeros(model.ndof,1); end
if nargin<6, overrides=struct(); end
settings=pae.options(overrides); time=time(:); r=model.ndof;
assert(numel(time)>=2 && all(isfinite(time)) && all(diff(time)>0), ...
    'pae:Time','Time must be finite and strictly increasing.');
q=zeros(r,numel(time)); v=q; a=q; loads=zeros(r,numel(time)); CL=NaN(numel(time),1);
iterations=zeros(numel(time),1); residualNorm=zeros(numel(time),1);
assert(numel(q0)==r && numel(v0)==r && all(isfinite([q0(:);v0(:)])), ...
    'pae:State','Invalid initial state.');
q(:,1)=q0(:); v(:,1)=v0(:);
M=model.structure.M; C=model.structure.C; external=model.structure.externalForce;
state=pae.evaluate(model,q(:,1),v(:,1),flow_at(time(1)));
a(:,1)=M\(state.Qaero+external-C*v(:,1)-state.internalForce);
loads(:,1)=state.Qaero; CL(1)=state.aero.CL;
beta=settings.newmarkBeta; gamma=settings.newmarkGamma;
for k=2:numel(time)
    dt=time(k)-time(k-1);
    qpred=q(:,k-1)+dt*v(:,k-1)+dt^2*(.5-beta)*a(:,k-1);
    vpred=v(:,k-1)+dt*(1-gamma)*a(:,k-1);
    next=qpred; flow=flow_at(time(k)); converged=false;
    for iteration=1:settings.maxIterations
        [residual,state,anew,vnew]=dynamic_residual(next);
        tolerance=settings.absoluteTolerance+settings.relativeTolerance* ...
            max([norm(state.Qaero+external),norm(state.internalForce),norm(M*anew)]);
        if norm(residual)<=tolerance, converged=true; break; end
        tangent=pae.jacobian(@only_residual,next,settings.differenceStep);
        assert(rcond(tangent)>1e-12,'pae:Transient','Singular transient tangent at t=%.6g.',time(k));
        correction=-tangent\residual;
        accepted=false;
        for search=0:settings.maxLineSearch-1
            candidate=next+2^(-search)*correction;
            if norm(only_residual(candidate))<(1-1e-4*2^(-search))*norm(residual)
                next=candidate; accepted=true; break;
            end
        end
        if ~accepted, break; end
    end
    [residual,state,anew,vnew]=dynamic_residual(next);
    tolerance=settings.absoluteTolerance+settings.relativeTolerance* ...
        max([norm(state.Qaero+external),norm(state.internalForce),norm(M*anew)]);
    converged=converged || norm(residual)<=tolerance;
    assert(converged,'pae:Transient','Transient equilibrium failed at t=%.6g, residual=%.6g.',time(k),norm(residual));
    q(:,k)=next; v(:,k)=vnew; a(:,k)=anew;
    loads(:,k)=state.Qaero; CL(k)=state.aero.CL;
    iterations(k)=iteration; residualNorm(k)=norm(residual);
end
response=struct('time',time,'q',q,'velocity',v,'acceleration',a,'Qaero',loads, ...
    'CL',CL,'iterations',iterations,'residualNorm',residualNorm, ...
    'structuralSource',model.structure.source,'fidelity','quasi-steady Newmark; no wake memory/added mass');
    function flow=flow_at(t)
        if isa(flowInput,'function_handle'), flow=flowInput(t); else, flow=flowInput; end
    end
    function [residual,evaluated,acceleration,velocity]=dynamic_residual(x)
        acceleration=(x-qpred)/(beta*dt^2);
        velocity=vpred+gamma*dt*acceleration;
        evaluated=pae.evaluate(model,x,velocity,flow);
        residual=M*acceleration+C*velocity+evaluated.internalForce-evaluated.Qaero-external;
    end
    function residual=only_residual(x)
        residual=dynamic_residual(x);
    end
end
