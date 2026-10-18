function result = external(count,bounds,Re,mask,tEnd,output,settings)
%EXTERNAL Laminar external-flow pilot using the own MAC/Brinkman core.
% Dimensionless U=c=rho=1. Periodic box plus explicitly declared fringe;
% not an unbounded-domain or turbulent DNS claim. Resolve space, eta and box.
if nargin<7,settings=struct();end
o=struct('referenceArea',1,'penaltyFactor',.05,'initialState',[],'saveEvery',100,'wallVelocity',[],'penaltyProjection','coupled','timeStepScale',1);
names=fieldnames(settings);for k=1:numel(names),assert(isfield(o,names{k}),'nsflow:Options','Unknown setting.');o.(names{k})=settings.(names{k});end
if ~isfolder(output),mkdir(output);end
g=nsflow.grid(count,bounds);d=g.dimension;nu=1/Re;
physical=struct('nu',nu,'rho',1,'mask',mask,'penaltyTime',o.penaltyFactor*min(g.h)^2/nu, ...
    'freestream',[1,zeros(1,d-1)],'penaltyProjection',o.penaltyProjection);
physical.fringe=@fringe;
if ~isempty(o.wallVelocity),physical.wallVelocity=o.wallVelocity;end
if isempty(o.initialState)
    U=cell(1,d);for i=1:d,U{i}=zeros(count);end
    U{1}=1-mask(g.velocityCoordinates(1,:),0);[U,p]=nsflow.project(g,U,.01,1);state=struct('U',{U},'p',p,'time',0);
else,state=o.initialState;end
initialTime=state.time;statisticsStart=initialTime+.6*(tEnd-initialTime);assert(tEnd>initialTime,'nsflow:Time','End time must exceed the initial state time.');
sumU=cell(1,d);sumU2=sumU;for i=1:d,sumU{i}=zeros(count);sumU2{i}=zeros(count);end
sumP=zeros(count);weight=0;trace=[];step=0;maxDivergence=0;maxMomentum=0;timer=tic;
while state.time<tEnd
    speed=0;for i=1:d,speed=speed+max(abs(state.U{i}),[],'all')/g.h(i);end
    dt=min([o.timeStepScale*.35/max(speed,1e-10),o.timeStepScale*.35/(nu*sum(1./g.h.^2)),tEnd-state.time]);
    [state,info]=nsflow.step(g,state,dt,physical);step=step+1;
    maxDivergence=max(maxDivergence,info.maximumDivergence);maxMomentum=max(maxMomentum,norm(info.momentumBalance));
    trace(end+1,:)=[state.time,dt,2*info.bodyForce/o.referenceArea,info.kineticEnergy,info.maximumDivergence]; %#ok<AGROW>
    if state.time>=statisticsStart
        for i=1:d,sumU{i}=sumU{i}+dt*state.U{i};sumU2{i}=sumU2{i}+dt*state.U{i}.^2;end
        sumP=sumP+dt*state.p;weight=weight+dt;
    end
    if mod(step,500)==0,fprintf('NS %dD %s: t*=%.3f, step %d, Cd=%.4f\n',d,mat2str(count),state.time,step,trace(end,3));end
end
meanU=sumU;fluctuationRMS=sumU;
for i=1:d,meanU{i}=sumU{i}/weight;fluctuationRMS{i}=sqrt(max(sumU2{i}/weight-meanU{i}.^2,0));end
meanP=sumP/weight;
result=struct('solidMask',mask(g.pressureCoordinates,state.time),'grid',g,'state',state,'meanVelocity',{meanU},'fluctuationRMS',{fluctuationRMS},'meanPressure',meanP,'trace',trace);
use=trace(:,1)>=statisticsStart;weights=trace(use,2);coeff=trace(use,3:2+d);average=sum(weights.*coeff,1)/sum(weights);
summary=struct('dimension',d,'count',count,'bounds',bounds,'Re',Re,'steps',step,'timeEnd',tEnd, ...
    'initialTime',initialTime,'statisticsWindow',[statisticsStart,tEnd],'meanForceCoefficient',average,'penaltyTime',physical.penaltyTime,'penaltyProjection',o.penaltyProjection,'timeStepScale',o.timeStepScale, ...
    'maximumDivergence',maxDivergence,'maximumMomentumDefect',maxMomentum,'seconds',toc(timer), ...
    'scope','own laminar Navier-Stokes; periodic MAC grid plus fringe and approximate Brinkman no-slip', ...
    'limitations','space/time/penalty/domain convergence required; no turbulence model and no tunnel-Re validation');
cad.write_json(fullfile(output,'summary.json'),summary);writematrix(trace,fullfile(output,'forces.csv'));
save(fullfile(output,'flow.mat'),'result','summary','-v7.3');result.summary=summary;
    function rate=fringe(coords)
        x=coords{1};lo=bounds(1,1);hi=bounds(1,2);width=.12*(hi-lo);
        distance=min(x-lo,hi-x);z=max(0,min(1,1-distance/width));rate=4*z.^2.*(3-2*z);
    end
end
