function results = run_tail_tests()
%RUN_TAIL_TESTS Independent plate, transfer and forced-dynamics verification.
names=["plate_rigid_modes","thin_cantilever","simply_supported_plate", ...
    "pressure_patch_conservation","harmonic_oscillator","forced_time_convergence", ...
    "pressure_history_bounds","modal_solid_adapter","published_data_audit","observable_comparison"];
checks={@rigid_modes,@thin_cantilever,@navier,@pressure_conservation,@harmonic, ...
    @forced_time,@pressure_history,@solid_adapter,@published_data,@observables};
passed=false(numel(names),1); metric=NaN(numel(names),1); details=strings(numel(names),1);
for k=1:numel(names)
    try
        [metric(k),details(k)]=checks{k}(); passed(k)=true;
    catch err
        details(k)=string(err.identifier)+": "+string(err.message);
    end
    fprintf('%s: %s | %s\n',names(k),string(passed(k)),details(k));
end
results=table(names',passed,metric,details,'VariableNames',{'test','passed','metric','details'});
assert(all(passed),'ael:Tests','Trailing-edge verification failed.');
end

function [metric,detail]=rigid_modes()
X=[0,0;.3,.02;.28,.2;-.01,.18];
[K,M]=plate.mitc4(X,.001,2.1e9,.3,1180);
q=reshape([.01+.2*X(:,1)-.1*X(:,2),.2*ones(4,1),-.1*ones(4,1)]',[],1);
metric=norm(K*q)/norm(K,'fro')/norm(q);
ev=eig(K); assert(sum(abs(ev)<1e-10*norm(K,2))==3 && metric<1e-13 && min(eig(M))>0);
detail="three rigid plate motions, no extra zero modes, positive consistent mass";
end

function [metric,detail]=thin_cantilever()
E=2.1e9; L=.2; errors=zeros(3,1); slenderness=[.01,.001,.0001];
for k=1:3
    h=L*slenderness(k); m=plate.mesh(linspace(0,L,33),linspace(0,.04,5),h);
    s=plate.assemble(m,struct('E',E,'nu',0,'rho',1180));
    fixed=fem3d.dofs(find(m.nodes(:,1)==0)); O=plate.probe(m,[L,.02],1);
    r=ael.static(s.K,fixed,s.P*ones(size(m.quads,1),1),O);
    exact=L^4/(8*(E*h^3/12)); errors(k)=abs(r.Y/exact-1);
end
metric=max(errors); assert(metric<.003 && errors(end)<.003);
detail=sprintf('uniform pressure beam limit nu=0, h/L=1e-2..1e-4, errors %.3g %.3g %.3g',errors);
end

function [metric,detail]=navier()
a=.2; b=.12; h=.001; E=2.1e9; nu=.3; rho=1180;
m=plate.mesh(linspace(0,a,21),linspace(0,b,13),h);
s=plate.assemble(m,struct('E',E,'nu',nu,'rho',rho)); X=m.nodes;
bx=find(X(:,1)==0 | X(:,1)==a); by=find(X(:,2)==0 | X(:,2)==b);
fixed=unique([3*[bx;by]-2;3*bx;3*by-1]);
model=ael.reduce(s.K,s.M,fixed,3,0,s.P,plate.probe(m,[a/2,b/2],1));
D=E*h^3/(12*(1-nu^2)); exact=pi/2*sqrt(D/(rho*h))*(1/a^2+1/b^2);
metric=abs(model.frequencyHz(1)/exact-1); assert(metric<.02);
detail="first simply-supported frequency against Navier closed form";
end

function [metric,detail]=pressure_conservation()
m=plate.mesh([0,.07,.13,.2],[0,.04,.1],.002);
xe=[0,.03,.16,.2]; ye=[0,.025,.08,.1]; B=plate.pressure_map(m,xe,ye);
p=(1:9)'; F=B*p; w=.7+.2*m.nodes(:,1)-.3*m.nodes(:,2);
[dx,dy]=ndgrid(diff(xe),diff(ye)); area=dx(:).*dy(:);
[cx,cy]=ndgrid((xe(1:end-1)+xe(2:end))/2,(ye(1:end-1)+ye(2:end))/2);
exact=[sum(area.*p);sum(area.*p.*cx(:));sum(area.*p.*cy(:))];
actual=[sum(F(1:3:end));m.nodes(:,1)'*F(1:3:end);m.nodes(:,2)'*F(1:3:end)];
work=abs(w'*F(1:3:end)-sum((.7+.2*cx(:)-.3*cy(:)).*area.*p));
metric=max([abs(actual-exact);work]); assert(metric<1e-13 && nnz(F(2:3:end))==0);
detail="nonmatching physical pressure patches conserve force, first moments and affine work";
end

function model=oscillator()
model=struct('M',2,'K',200,'C',.8,'B',1,'O',1);
end

function [metric,detail]=harmonic()
model=oscillator(); f=[0,1,10/(2*pi),3]; r=ael.harmonic(model,f,3);
exact=3./(200-2*(2*pi*f).^2+1i*.8*2*pi*f);
metric=max(abs(r.Y-exact)); assert(metric<1e-13);
detail="complex amplitude and phase against exact damped oscillator at resonance and DC";
end

function [metric,detail]=forced_time()
model=oscillator(); w=8; H=1/(200-2*w^2+1i*.8*w); errors=zeros(2,1);
for k=1:2
    t=linspace(0,3,300*2^k+1);
    r=ael.transient(model,t,@(t) cos(w*t),real(H),real(1i*w*H));
    exact=real(H*exp(1i*w*t)); errors(k)=max(abs(r.Y-exact))/abs(H);
end
metric=errors(2); assert(metric<.002 && errors(2)<.3*errors(1));
detail=sprintf('exact forced steady state, second-order time convergence %.3g -> %.3g',errors);
end

function [metric,detail]=pressure_history()
fun=ael.pressure_history([0,1,2],[1,3,5;4,2,0]);
metric=norm(fun(.5)-[2;3]); caught=false;
try, fun(2.01); catch err, caught=strcmp(err.identifier,'ael:PressureTime'); end
assert(metric<1e-14 && caught); detail="pressure channel interpolation and extrapolation rejection";
end

function [metric,detail]=solid_adapter()
m=fem3d.box_mesh([.1,.02,.01],[8,2,2]); s=fem3d.assemble(m,fem3d.material_reference());
fixed=fem3d.dofs(find(m.nodes(:,1)==0)); original=fem3d.modal(s,struct('dofs',fixed),3);
B=sparse(size(s.K,1),1); B(end)=1; O=B'; model=ael.reduce(s.K,s.M,fixed,3,.01,B,O);
metric=max(abs(model.frequencyHz./original.frequencyHz-1));
assert(metric<1e-9 && norm(model.M-eye(3),'fro')<1e-8);
detail="same structural adapter accepts solid interleaved translations without plate assumptions";
end

function [metric,detail]=published_data()
cases=ael.reference_cases(); assert(numel(cases)==9);
count=0;
for k=1:numel(cases)
    T=cases(k).reference; assert(height(T)==3 && all(T.omega_nd_mean>0) && ...
        all(T.omega_nd_sample_sd>0) && T.source_mode_label(1)=="1-B");
    count=count+sum(T.published_mean_uses_all_three==0);
end
metric=count; assert(count==6);
detail="nine geometry keys, 27 measured observations; six published subset means exposed";
end

function [metric,detail]=observables()
t=(0:999)/100; Y=2+.3*sin(2*pi*5*t); referenceT=t(1:2:end);
r=ael.compare_series(t,Y,referenceT,Y(1:2:end)+.01,.25);
m=ael.response_metrics(t,Y);
metric=abs(r.rmse-.01)+abs(m.dominantFrequencyHz-5);
assert(metric<1e-12 && abs(m.fluctuationRMS-.3/sqrt(2))<1e-12 && ...
    abs(r.normalizedRMSE-.04)<1e-12);
detail="known bias/RMSE and prescribed 5 Hz displacement, physical normalization and RMS";
end
