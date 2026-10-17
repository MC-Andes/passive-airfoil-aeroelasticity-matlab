function report = run_viscous_tests()
%RUN_VISCOUS_TESTS Own 2D/3D incompressible MAC operators against exact solutions.
names=["mac_projection","poiseuille_channel","taylor_green_decay","abc_three_dimensional","immersed_pressure_momentum"];
checks={@projection,@channel,@taylor,@abc,@immersed};metric=zeros(numel(names),1);details=strings(numel(names),1);
for k=1:numel(names),[metric(k),details(k)]=checks{k}();fprintf('%s: true | %s\n',names(k),details(k));end
report=table(names',true(numel(names),1),metric,details,'VariableNames',{'test','passed','metric','details'});
end
function [v,d]=projection()
g=nsflow.grid([16,20],[0,2*pi;0,2*pi]);rng(91);U={randn(16,20),randn(16,20)};[~,~,v]=nsflow.project(g,U,.01,1);assert(v<1e-12);d="periodic pressure solve eliminates discrete mass divergence to roundoff";
end
function [v,d]=channel()
e=zeros(3,1);for k=1:3,r=nsflow.channel(8*2^(k-1),1,.1,1,1);e(k)=r.relativeError;end
v=e(end);assert(all(e(1:end-1)./e(2:end)>3.8) && v<.001);d=sprintf('no-slip Poiseuille errors %.3g %.3g %.3g',e);
end
function [v,d]=taylor()
e=zeros(3,1);p=struct('nu',.1,'rho',1);
for k=1:3
 n=8*2^(k-1);g=nsflow.grid([n,n],[0,2*pi;0,2*pi]);[U,P]=nsflow.analytic(g,0,p.nu,'taylor_green');s=struct('U',{U},'p',P,'time',0);dt=.0025;
 for j=1:40,[s,info]=nsflow.step(g,s,dt,p);end
 [exact,~]=nsflow.analytic(g,s.time,p.nu,'taylor_green');e(k)=sqrt(sum((s.U{1}-exact{1}).^2+(s.U{2}-exact{2}).^2,'all')/sum(exact{1}.^2+exact{2}.^2,'all'));
 assert(info.maximumDivergence<1e-12);
end
v=e(end);assert(e(1)/e(2)>3.5 && e(2)/e(3)>3.5 && v<1e-3);d=sprintf('exact 2D viscous vortex decay, spatial errors %.3g %.3g %.3g',e);
end
function [v,d]=abc()
e=zeros(2,1);p=struct('nu',.1,'rho',1);
for k=1:2
 n=8*2^(k-1);g=nsflow.grid([n,n,n],repmat([0,2*pi],3,1));[U,P]=nsflow.analytic(g,0,p.nu,'abc');s=struct('U',{U},'p',P,'time',0);
 for j=1:20,[s,info]=nsflow.step(g,s,.0025,p);end
 exact=nsflow.analytic(g,s.time,p.nu,'abc');error=0;scale=0;
 for i=1:3,error=error+sum((s.U{i}-exact{i}).^2,'all');scale=scale+sum(exact{i}.^2,'all');end,e(k)=sqrt(error/scale);assert(info.maximumDivergence<1e-12);
end
v=e(end);assert(v<.002 && e(1)/e(2)>3.5);d=sprintf('exact 3D Beltrami/ABC Navier-Stokes decay, errors %.3g %.3g',e);
end

function [v,d]=immersed()
g=nsflow.grid([32,32],[-2,2;-2,2]);rng(71);U={randn(32),randn(32)};p=struct('rho',1,'nu',.05,'mask',@(coords,t)nsflow.circle_mask(coords,.04),'penaltyTime',.001);[V,P,info]=nsflow.project_penalized(g,U,0,.01,p);assert(info.maximumDivergence<1e-7);s=struct('U',{V},'p',P,'time',0);[~,check]=nsflow.step(g,s,.0001,p);v=norm(check.momentumBalance);assert(v<1e-10);
g2=nsflow.grid([24,24],[-pi,pi;-pi,pi]);X=g2.velocityCoordinates;R={.01*2*sin(g2.h(1)/2)/g2.h(1)*cos(X{1,1}).*cos(X{1,2}),-.01*2*sin(g2.h(2)/2)/g2.h(2)*sin(X{2,1}).*sin(X{2,2})};[quiet,~,~]=nsflow.project_penalized(g2,R,0,.01,p);assert(max(abs(quiet{1}),[],'all')<1e-8 && max(abs(quiet{2}),[],'all')<1e-8);d="coupled Brinkman/pressure projection, mean-zero pressure and exact body-reaction/global momentum balance";
end
