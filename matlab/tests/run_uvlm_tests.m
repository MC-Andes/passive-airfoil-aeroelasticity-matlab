function report = run_uvlm_tests()
%RUN_UVLM_TESTS Independent rigid/unsteady aerodynamics and coupling checks.
names=["biot_savart_ring","uvlm_steady_limit","kelvin_and_work","harmonic_theodorsen","free_wake","static_feedback","dynamic_energy","nonlinear_tangent_consistency","wake_vertex_derivative","wake_time_subdivision"];
checks={@kernel,@steady,@kelvin_work,@harmonic,@free_wake,@static_feedback,@dynamic_energy,@tangent_consistency,@wake_derivative,@time_subdivision};
metric=zeros(numel(names),1);details=strings(numel(names),1);
for k=1:numel(names),[metric(k),details(k)]=checks{k}();fprintf('%s: true | %s\n',names(k),details(k));end
report=table(names',true(numel(names),1),metric,details,'VariableNames',{'test','passed','metric','details'});
end
function [v,d]=time_subdivision()
[s,f]=flat(4,3,3);s.nodes(:,:,3)=.02*sin(pi*s.nodes(:,:,1));phi=zeros(numel(s.nodes),2);phi(3:3:end,1)=1;phi(3:3:end,2)=reshape(s.nodes(:,:,1),[],1);
a=uvlm.linear_model(s.nodes,phi,f,.02,51,struct('wakeLength',1000,'compactTangent',true));b=uvlm.refine_time_model(a);c=uvlm.linear_model(s.nodes,phi,f,.01,101,struct('wakeLength',1000,'compactTangent',true));
fields={'W','Fw','Gq','Fq','D','projectedD','Lq'};errors=zeros(1,numel(fields));for k=1:numel(fields),x=b.(fields{k});y=c.(fields{k});errors(k)=norm(x(:)-y(:))/max(norm(y(:)),1e-15);end
assert(max(errors(1:4))<1e-9 && max(errors(5:end))<1e-6);m=struct('K',diag([100,150]),'M',diag([10,2]),'C',eye(2)*.1);time=0:.01:1;r=fsi.linear_dynamic(m,b,time,[1e-5;0],[0;0]);z=fsi.linear_dynamic(m,c,time,[1e-5;0],[0;0]);v=norm(r.Q-z.Q,'fro')/norm(z.Q,'fro');assert(v<1e-6);d="exact half-step ring/vertex subdivision and first-row shape correction against independent direct kernels and full coupled trajectories on a cambered lattice";
end
function [v,d]=wake_derivative()
rng(81);P=randn(30,3)*.1;segment=[.001,0,.00005];centers=[.2,.1,.01]+[1;10;100;1000;10000]*segment;
D=uvlm.wake_vertex_gradient(P,centers,segment,1e-5);errors=zeros(1,3);
for axis=1:3
 h=1e-20;cp=complex(centers);cp(:,axis)=cp(:,axis)+1i*h;
 exact=imag(uvlm.segment_influence(P,centers-segment,cp,1e-5)+uvlm.segment_influence(P,cp,centers+segment,1e-5))/h;
 errors(axis)=norm(D{axis}(:)-exact(:))/norm(exact(:));
end
v=max(errors);assert(v<1e-9);d="analytic shared wake vertex derivative against independent complex-step finite-segment kernels over five age scales";
end
function [s,f]=flat(nx,ny,aspect)
p=struct('spanStations',[-aspect/2;aspect/2],'leadingEdge',[0,-aspect/2,0;0,aspect/2,0],'chord',[1;1]);s=pae.make_surface(p,nx,ny);s.nodes(:,:,1)=repmat(linspace(0,1,nx+1)',1,ny+1);
f=struct('rho',1,'speed',1,'alpha',.03,'referenceArea',aspect,'referenceChord',1);
end
function [v,d]=kernel()
vertices=[0,0,0;0,1,0;1,1,0;1,0,0];P=[.2,.4,.1;.7,.3,-.2];expected=zeros(size(P));
for k=1:4,expected=expected+pae.segment_velocity(P,vertices(k,:),vertices(mod(k,4)+1,:),1e-5);end
v=norm(uvlm.ring_velocity(P,vertices,1:4,1,1e-5)-expected,'fro');assert(v<1e-13);d="vectorized ring kernel against four independent finite-segment evaluations";
end
function [v,d]=steady()
[s,f]=flat(6,10,6);[state,r]=uvlm.initialize(s.nodes,zeros(numel(s.nodes)/3,3),f,struct('wakeLength',100));a=pae.vlm(s.nodes,zeros(numel(s.nodes)/3,3),f);
v=abs(r.CL/a.CL-1);assert(v<.01);
[state,next]=uvlm.step(state,s.nodes,zeros(numel(s.nodes)/3,3),f,.05);assert(abs(next.CL/r.CL-1)<1e-4 && next.boundaryResidual<1e-12);
d="stationary ring/wake limit against independent horseshoe VLM, including persistence in time";
end
function [v,d]=kelvin_work()
[s,f]=flat(3,3,3);g=uvlm.lattice(s.nodes);rng(81);gamma=randn(9,1)*.01;wake=struct('nodes',cat(1,reshape(g.trailingEdge,1,4,3),reshape(g.trailingEdge+[2,0,0],1,4,3)),'gamma',zeros(1,3));
r=uvlm.loads(g,gamma,gamma*.9,.1,wake,f,0,zeros(numel(s.nodes)/3,3),uvlm.options());
[state,~]=uvlm.initialize(s.nodes,zeros(numel(s.nodes)/3,3),f);[next,out]=uvlm.step(state,s.nodes,zeros(numel(s.nodes)/3,3),f,.02);
assert(out.kelvinDefect==0 && isequal(next.wake.gamma(2:end,:),state.wake.gamma));
% Verify global force/moment under an independent affine velocity field.
X=reshape(s.nodes,[],3);omega=[.1,.2,-.3];velocity=[.2,-.1,.4]+cross(repmat(omega,size(X,1),1),X,2);
work=sum(velocity.*r.nodalForces,'all');expected=dot([.2,-.1,.4],r.totalForce)+dot(omega,r.totalMoment);
v=abs(work-expected);assert(v<1e-12);d="immutable old wake, Kutta shedding and affine force/moment/virtual-work balance";
end
function [v,d]=harmonic()
% Modest verification mesh; the campaign refines space/time further.
[s,f]=flat(16,1,1000);f.alpha=0;k=.5;h=.001;T=2*pi;dt=1/64;time=0:dt:4*T;
velocity=zeros(numel(s.nodes)/3,3);X=s.nodes;X(:,:,3)=h;[state,~]=uvlm.initialize(X,velocity,f,struct('wakeLength',100,'maximumWakeRows',8000));cl=zeros(size(time));
for j=2:numel(time),X=s.nodes;X(:,:,3)=h*cos(time(j));velocity(:,3)=-h*sin(time(j));[state,r]=uvlm.step(state,X,velocity,f,dt);cl(j)=r.CL;end
use=time>2*T;fit=[cos(time(use))',sin(time(use))',ones(sum(use),1)]\cl(use)';observed=fit(1)-1i*fit(2);
C=besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k));expected=2*pi*k^2*h-4*pi*1i*k*C*h;
v=abs(observed/expected-1);phase=abs(rad2deg(angle(observed/expected)));assert(v<.04 && phase<2);
d=sprintf('Theodorsen heave k=0.5: complex error %.3g, phase %.3g deg; tighter campaign required for closure',v,phase);
end
function [v,d]=free_wake()
[s,f]=flat(2,3,3);[state,~]=uvlm.initialize(s.nodes,zeros(numel(s.nodes)/3,3),f,struct('convection','free','wakeLength',5));
for k=1:4,[state,r]=uvlm.step(state,s.nodes,zeros(numel(s.nodes)/3,3),f,.05);end
v=max(abs(state.wake.nodes(2:end,:,3)),[],'all');assert(v>1e-7 && all(isfinite(state.wake.nodes(:))) && r.kelvinDefect==0);d="wake vertices convect with computed induced velocity, rather than decorative streamlines";
end
function [v,d]=static_feedback()
X=zeros(2,2,3);model=struct('X0',X,'PhiSurface',zeros(12,1),'K',10);model.PhiSurface(3:3:end)=.5;model.aero=@force;
r=fsi.static(model,struct(),struct('relaxation',.6,'loadTolerance',1e-8,'displacementTolerance',1e-8));expected=1/(10-.2);v=abs(r.q/expected-1);assert(v<1e-7);d="static two-way feedback against an independently solved scalar aeroelastic equilibrium";
 function a=force(nodes,velocity,flow) %#ok<INUSD>
  displacement=mean(nodes(:,:,3),'all');load=1+.4*displacement;a=struct('nodalForces',repmat([0,0,load/2],4,1));
 end
end
function [v,d]=dynamic_energy()
model=struct('X0',zeros(2,2,3),'PhiSurface',zeros(12,1),'K',100,'M',1,'C',.1);model.PhiSurface(3:3:end)=.5;
model.initializeFluid=@initialize;model.advanceFluid=@advance;time=linspace(0,1,401);r=fsi.dynamic(model,struct(),time,.001,0,struct('relativeTolerance',1e-8));
v=r.relativeEnergyBalance;assert(r.converged && v<1e-7 && r.maximumInterfaceWorkDefect<1e-14);d="strong partitioned Newmark energy and interface work under independently specified zero fluid forcing";
 function [state,aero]=initialize(nodes,velocity,flow) %#ok<INUSD>
 state=struct('time',0);aero=struct('nodalForces',zeros(4,3),'CL',0);
 end
 function [state,aero]=advance(old,nodes,velocity,flow,dt) %#ok<INUSD>
 state=struct('time',old.time+dt);aero=struct('nodalForces',zeros(4,3),'CL',0);
 end
end

function [v,d]=tangent_consistency()
p=struct('spanStations',[-1.5;1.5],'leadingEdge',[0,-1.5,0;0,1.5,0],'chord',[1;1]);s=pae.make_surface(p,4,3);s.nodes(:,:,1)=repmat(linspace(0,1,5)',1,4);f=struct('rho',1,'speed',1,'alpha',.03,'referenceChord',1,'referenceArea',3);phi=zeros(numel(s.nodes),1);phi(3:3:end)=1;dt=.02;time=0:dt:2;fluid=uvlm.linear_model(s.nodes,phi,f,dt,numel(time),struct('wakeLength',1000));[state,base]=uvlm.initialize(s.nodes,zeros(numel(s.nodes)/3,3),f,struct('wakeLength',1000));gam=zeros(12,1);gh=zeros(fluid.spanPanels*fluid.wakeRows,1);qh=zeros(fluid.wakeRows,1);lastq=0;out=[];
for j=2:numel(time)
 gh=[fluid.TEextract*gam;gh(1:end-fluid.spanPanels)];qh=[lastq;qh(1:end-1)];q=1e-5*(1-cos(4*time(j)));v=4e-5*sin(4*time(j));gnew=fluid.D*q+fluid.E*v+fluid.A\(-fluid.W*gh-fluid.Gq*qh);
 force=fluid.load0+fluid.Lg*gnew+fluid.Ld/dt*(gnew-gam)+fluid.Lq*q+fluid.Lv*v+fluid.Fw*gh+fluid.Fq*qh;
 X=s.nodes;X(:,:,3)=q;vel=zeros(numel(s.nodes)/3,3);vel(:,3)=v;[state,r]=uvlm.step(state,X,vel,f,dt);full=phi'*reshape(r.nodalForces',[],1);out(end+1,:)=[time(j),force,full];gam=gnew;lastq=q;
end
error=max(abs(out(:,2)-out(:,3)))/max(abs(out(:,3)-fluid.load0));v=error;assert(v<1e-4);direction=sin((1:numel(fluid.gamma0))');h=1e-5;gp=fluid.gamma0+h*direction;gm=fluid.gamma0-h*direction;ap=uvlm.loads(fluid.state0.geometry,gp,gp,Inf,fluid.state0.wake,f,0,zeros(numel(s.nodes)/3,3),fluid.settings);am=uvlm.loads(fluid.state0.geometry,gm,gm,Inf,fluid.state0.wake,f,0,zeros(numel(s.nodes)/3,3),fluid.settings);fd=phi'*(reshape(ap.nodalForces',[],1)-reshape(am.nodalForces',[],1))/(2*h);assert(norm(fd-fluid.Lg*direction)/max(norm(fd),1e-12)<1e-9);compact=uvlm.linear_model(s.nodes,phi,f,dt,numel(time),struct('wakeLength',1000,'compactTangent',true));m=struct('K',100,'M',10,'C',.1);a=fsi.linear_dynamic(m,fluid,time,1e-5,0,1);b=fsi.linear_dynamic(m,fluid,time,1e-5,0,16);c=fsi.linear_dynamic(m,compact,time,1e-5,0,16);assert(norm(c.Q-b.Q,'fro')/norm(b.Q,'fro')<1e-11);assert(norm(a.Q-b.Q,'fro')/norm(a.Q,'fro')<1e-11);d="nonlinear time-marching UVLM versus its independently assembled tangent for a prescribed small motion";

end
