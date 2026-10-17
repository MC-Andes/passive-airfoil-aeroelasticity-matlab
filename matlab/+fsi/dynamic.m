function result = dynamic(model,flow,time,q0,v0,settings)
%DYNAMIC Strong partitioned Newmark/UVLM with material work-conjugate exchange.
% Every subiteration starts from the same OLD fluid state. Wake is committed
% only once per converged time step. No prescribed oscillating aerodynamic load.
if nargin<6,settings=struct();end
o=struct('maximumIterations',50,'relaxation',.7,'relativeTolerance',1e-6, ...
    'forceScale',1e-8,'errorOnFailure',true,'uvlmOptions',struct(),'storePressureHistory',false,'maximumPressureBytes',268435456);
names=fieldnames(settings); for k=1:numel(names),assert(isfield(o,names{k}),'fsi:Options','Unknown setting.');o.(names{k})=settings.(names{k});end
time=time(:)';nt=numel(time); n=size(model.K,1);
assert(nt>=2 && all(diff(time)>0) && time(1)==0,'fsi:Time','Time must increase from zero.');
Q=zeros(n,nt);V=Q;A=Q;F=Q;Q(:,1)=q0;V(:,1)=v0;
X=geometry(q0);velocity=reshape(model.PhiSurface*v0,3,[])';
if isfield(model,'initializeFluid'),[fluid,aero]=model.initializeFluid(X,velocity,flow);
else,[fluid,aero]=uvlm.initialize(X,velocity,flow,o.uvlmOptions);end
Fa=reshape(aero.nodalForces',[],1);F(:,1)=model.PhiSurface'*Fa;
if isfield(model,'nodalAccelerationLoad')
    added=model.PhiSurface'*model.nodalAccelerationLoad;A(:,1)=(model.M-added)\(F(:,1)-model.C*v0-model.K*q0);initialAddedForce=model.nodalAccelerationLoad*A(:,1);Fa=Fa+initialAddedForce;F(:,1)=model.PhiSurface'*Fa;
else
    A(:,1)=model.M\(F(:,1)-model.C*v0-model.K*q0);
end
pressure=[];if o.storePressureHistory
 assert(isfield(aero,'pressureJump') && numel(aero.pressureJump)*nt*8<=o.maximumPressureBytes,'fsi:Pressure','Pressure history is unavailable or exceeds the explicit memory budget.');
 pressure=zeros(numel(aero.pressureJump),nt);pressure(:,1)=aero.pressureJump;
 if isfield(model,'pressureAccelerationLoad'),pressure(:,1)=pressure(:,1)+model.pressureAccelerationLoad*A(:,1);end
end
energy=zeros(1,nt);work=energy;dissipation=energy;balance=energy;exchange=energy;
energy(1)=.5*(v0'*model.M*v0+q0'*model.K*q0);
iterations=zeros(1,nt);residual=zeros(1,nt);CL=zeros(1,nt);CL(1)=aero.CL;
if isfield(model,'nodalAccelerationLoad') && all(isfield(flow,{'rho','speed','alpha','referenceArea'})) && flow.speed>0
 referenceFlow=flow;if isfield(referenceFlow,'gust'),referenceFlow=rmfield(referenceFlow,'gust');end
 direction=uvlm.air_velocity(referenceFlow,0,[0,0,0]);direction=direction/norm(direction);lift=cross(direction,[0,1,0]);lift=lift/norm(lift);
 CL(1)=CL(1)+dot(sum(reshape(initialAddedForce,3,[])',1),lift)/(.5*flow.rho*flow.speed^2*flow.referenceArea);
end
snapshots={};snapshotTime=[];oldDt=NaN;
for step=2:nt
    dt=time(step)-time(step-1); up=Q(:,step-1)+dt*V(:,step-1)+.25*dt^2*A(:,step-1);
    vp=V(:,step-1)+.5*dt*A(:,step-1); guess=A(:,step-1);converged=false;
    if dt~=oldDt,effective=model.M+.5*dt*model.C+.25*dt^2*model.K;oldDt=dt;end
    for iteration=1:o.maximumIterations
        q=up+.25*dt^2*guess;v=vp+.5*dt*guess;
        X=geometry(q); velocity=reshape(model.PhiSurface*v,3,[])';
        if isfield(model,'advanceFluid'),[trial,aero]=model.advanceFluid(fluid,X,velocity,flow,dt);
        else,[trial,aero]=uvlm.step(fluid,X,velocity,flow,dt);end
        nextFa=reshape(aero.nodalForces',[],1);load=model.PhiSurface'*nextFa;
        defect=model.M*guess+model.C*v+model.K*q-load;
        relative=norm(defect)/max([norm(load),norm(model.K*q),norm(model.M*guess),o.forceScale]);
        if relative<=o.relativeTolerance,converged=true;break;end
        target=effective\(load-model.C*vp-model.K*up);
        guess=guess+o.relaxation*(target-guess);
    end
    if ~converged
        if o.errorOnFailure,error('fsi:Dynamic','Subiterations failed at t=%g s; residual %g.',time(step),relative);end
        % Never return an unconverged trace as a successful simulation.
        result=struct('converged',false,'failedStep',step,'time',time(1:step-1),'relativeResidual',relative);return;
    end
    Q(:,step)=q;V(:,step)=v;A(:,step)=guess;F(:,step)=load;
    if o.storePressureHistory,pressure(:,step)=aero.pressureJump;end
    dq=q-Q(:,step-1);vmid=(v+V(:,step-1))/2;
    increment=.5*(load+F(:,step-1))'*dq;
    exchange(step)=abs(increment-.5*(Fa+nextFa)'*(model.PhiSurface*dq));
    work(step)=work(step-1)+increment;dissipation(step)=dissipation(step-1)+dt*vmid'*model.C*vmid;
    energy(step)=.5*(v'*model.M*v+q'*model.K*q);
    balance(step)=energy(step)-energy(1)-work(step)+dissipation(step);
    iterations(step)=iteration;residual(step)=relative;CL(step)=aero.CL;fluid=trial;Fa=nextFa;
    if isfield(fluid,'wake') && (mod(step,max(1,floor(nt/10)))==0 || step==nt)
        snapshots{end+1}=fluid.wake;snapshotTime(end+1)=time(step); %#ok<AGROW>
    end
end
scale=max([max(abs(energy)),max(abs(work)),max(abs(dissipation)),1e-15]);
result=struct('converged',true,'time',time,'Q',Q,'V',V,'A',A,'force',F,'CL',CL, ...
    'energy',energy,'aerodynamicWork',work,'dissipation',dissipation,'energyBalance',balance, ...
    'relativeEnergyBalance',max(abs(balance))/scale,'maximumInterfaceWorkDefect',max(exchange), ...
    'iterations',iterations,'relativeResidual',residual,'fluid',fluid, ...
    'wakeSnapshots',{snapshots},'wakeTimes',snapshotTime,'pressureHistory_Pa',pressure,'settings',o);
    function X=geometry(q)
        x=reshape(model.X0,[],3)+reshape(model.PhiSurface*q,3,[])';X=reshape(x,size(model.X0));
    end
end
