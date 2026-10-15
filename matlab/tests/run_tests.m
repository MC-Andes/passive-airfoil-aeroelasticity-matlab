function results = run_tests()
%RUN_TESTS Physics/reference checks, executed in MATLAB without extra toolboxes.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
checks={@rigid_symmetry,@thin_airfoil_limit,@coefficient_scaling,@mesh_convergence, ...
    @virtual_work,@analytical_feedback,@divergence_failure,@free_oscillator, ...
    @steady_transient,@affine_modal_interpolation,@reference_geometry};
names=string(cellfun(@func2str,checks,'UniformOutput',false))'; elapsed=zeros(numel(checks),1);
for i=1:numel(checks)
    started=tic; checks{i}(); elapsed(i)=toc(started);
    fprintf('PASS %-30s %.3f s\n',names(i),elapsed(i));
end
results=table(names,true(numel(checks),1),elapsed,'VariableNames',{'check','passed','seconds'});
end

function [s,f]=flat_surface(nx,ny,aspect)
p=struct('spanStations',[-aspect*.125;aspect*.125], ...
    'leadingEdge',[0,-aspect*.125,0;0,aspect*.125,0],'chord',[.25;.25]);
s=pae.make_surface(p,nx,ny);
f=struct('rho',1.225,'speed',15,'alpha',deg2rad(3), ...
    'referenceArea',s.referenceArea,'referenceChord',s.referenceChord, ...
    'referencePoint',[.25/4,0,0]);
end

function rigid_symmetry()
[s,f]=flat_surface(5,20,6); velocity=zeros(numel(s.nodes)/3,3);
positive=pae.vlm(s.nodes,velocity,f); f.alpha=-f.alpha;
negative=pae.vlm(s.nodes,velocity,f); f.alpha=0;
zero=pae.vlm(s.nodes,velocity,f);
assert(abs(positive.CL+negative.CL)<1e-10 && positive.CL>0);
assert(norm(zero.totalForce)<1e-10 && positive.CDi>0);
assert(positive.boundaryResidual<1e-10);
end

function thin_airfoil_limit()
[s,f]=flat_surface(8,48,100); f.alpha=deg2rad(1);
a=pae.vlm(s.nodes,zeros(numel(s.nodes)/3,3),f);
ratio=a.CL/f.alpha/(2*pi);
fprintf('  Finite AR=100 slope / (2*pi) = %.6f\n',ratio);
assert(ratio>.95 && ratio<1.04,'High-aspect-ratio lift slope inconsistent with thin-airfoil limit.');
end

function coefficient_scaling()
[s,f]=flat_surface(4,16,5); v=zeros(numel(s.nodes)/3,3);
a=pae.vlm(s.nodes,v,f); f.speed=2*f.speed; f.rho=.7*f.rho;
b=pae.vlm(s.nodes,v,f);
assert(abs(a.CL-b.CL)<1e-11 && abs(a.CDi-b.CDi)<1e-11);
assert(norm(b.totalForce-2.8*a.totalForce)<1e-10);
end

function mesh_convergence()
[s,f]=flat_surface(4,16,6); a=pae.vlm(s.nodes,zeros(numel(s.nodes)/3,3),f);
[s,f]=flat_surface(8,32,6); b=pae.vlm(s.nodes,zeros(numel(s.nodes)/3,3),f);
difference=abs(a.CL/b.CL-1); fprintf('  Lattice refinement lift difference = %.4f %%\n',100*difference);
assert(difference<.05,'Lift has not converged within 5%% for this reference refinement.');
end

function virtual_work()
[s,~]=flat_surface(2,4,3); N=numel(s.nodes)/3; X=reshape(s.nodes,[],3);
st=struct('M',1,'C',.1,'K',10,'Phi',zeros(3*N,1),'labels',{{'rotation'}},'source','test');
st.geometryFcn=@(q) rotated_geometry(X,q);
model=pae.make_model(s,st); q=.08; [nodes,J]=pae.geometry(model,q); p=pae.panels(nodes);
forces=[sin((1:numel(p.area))'),cos((1:numel(p.area))'),ones(numel(p.area),1)];
Q=pae.map_loads(J,p.Hforce,forces); h=1e-7;
[plus,~]=pae.geometry(model,q+h); [minus,~]=pae.geometry(model,q-h);
motion=p.Hforce*(reshape(plus,[],3)-reshape(minus,[],3))/(2*h);
work=sum(forces.*motion,'all'); assert(abs(Q-work)<1e-8*max(1,abs(work)));
end

function [X,J]=rotated_geometry(X0,q)
R=[cos(q),0,sin(q);0,1,0;-sin(q),0,cos(q)];
dR=[-sin(q),0,cos(q);0,0,0;-cos(q),0,-sin(q)];
X=X0*R'; J=reshape(X0*dR',[],1);
end

function model=feedback_model()
[s,~]=flat_surface(2,4,3); N=numel(s.nodes)/3;
basis=zeros(N,3); basis(:,3)=1;
st=struct('M',2,'C',.4,'K',100,'Phi',basis(:),'labels',{{'heave'}},'source','analytical reference');
st.aeroFcn=@analytical_aero; model=pae.make_model(s,st);
end

function a=analytical_aero(nodes,velocity,f)
p=pae.panels(nodes); height=mean(nodes(:,:,3),'all'); speed=mean(velocity(:,3));
load=f.force0+f.feedback*f.speed^2*height-f.aeroDamping*speed;
forces=zeros(numel(p.area),3); forces(:,3)=load/numel(p.area);
a=struct('panelForces',forces,'CL',NaN,'CDi',NaN,'Cm',NaN,'fidelity','analytical feedback test');
end

function f=feedback_flow()
f=struct('rho',1,'speed',10,'alpha',0,'force0',10,'feedback',.3,'aeroDamping',.8);
end

function analytical_feedback()
model=feedback_model(); f=feedback_flow(); solution=pae.solve_static(model,f);
expected=f.force0/(model.structure.K-f.feedback*f.speed^2);
assert(abs(solution.q-expected)<1e-10 && solution.residualNorm<1e-9);
local=pae.linearize(model,solution);
expectedA=[0,1;-(100-.3*100)/2,-(.4+.8)/2];
assert(norm(local.A-expectedA)<1e-8);
end

function divergence_failure()
model=feedback_model(); f=feedback_flow(); f.speed=sqrt(100/.3);
options=struct('errorOnFailure',false);
solution=pae.solve_static(model,f,[],options);
assert(~solution.converged && solution.residualNorm>1, ...
    'A singular static problem must not produce a reported equilibrium.');
end

function free_oscillator()
model=feedback_model(); f=feedback_flow(); f.speed=0; f.force0=0; f.feedback=0; f.aeroDamping=0;
time=(0:.002:2)'; q0=.01;
response=pae.simulate(model,f,time,q0,0);
M=2; K=100; C=.4; decay=C/(2*M); omega=sqrt(K/M-decay^2);
exact=q0*exp(-decay*time).*(cos(omega*time)+decay/omega*sin(omega*time));
error=max(abs(response.q'-exact)); fprintf('  Damped oscillator max absolute error = %.3g m\n',error);
assert(error<3e-6,'Newmark response disagrees with analytical damped oscillator.');
energy=.5*M*response.velocity.^2+.5*K*response.q.^2;
assert(all(diff(energy)<1e-11) && energy(end)<energy(1));
end

function steady_transient()
[model,f]=pae.demo_case(3,8); solution=pae.solve_static(model,f);
response=pae.simulate(model,f,(0:.01:.1)',solution.q,[]);
assert(max(abs(response.q-solution.q),[],'all')<1e-9);
assert(max(response.residualNorm)<1e-7);
end

function affine_modal_interpolation()
points=[0,0,0;1,0,0;0,1,0;1,1,0;.2,.4,0;.8,.3,0];
targets=[.1,.2,0;.5,.5,0;.9,.8,0];
modes=zeros(size(points)); modes(:,3)=1+.2*points(:,1)+.3*points(:,2);
Phi=pae.interpolate_basis(points,modes,targets);
field=reshape(Phi,[],3);
assert(max(abs(field(:,3)-(1+.2*targets(:,1)+.3*targets(:,2))))<1e-10);
end

function reference_geometry()
projectRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));
[surface,info]=pae.surface_from_sections(fullfile(projectRoot,'data','reference'),4,10,"constant");
assert(info.solid_count==1 && strcmp(info.units,'m'));
assert(abs(surface.span-.36)<1e-5 && abs(surface.referenceChord-.247)<.002);
assert(surface.referenceArea>.08 && surface.referenceArea<.10);
assert(all(isfinite(surface.nodes(:))));
end
