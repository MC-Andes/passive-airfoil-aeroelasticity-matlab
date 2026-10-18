function report = run_cad_aeroelastic(bundlePath,settings)
%RUN_CAD_AEROELASTIC Conditional CAD one/two-way statics and UVLM dynamics.
% Uses actual CAD-derived structural modes/transfer, never a replacement plate.
if nargin<2,settings=struct();end
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'matlab'));
o=struct('duration',2,'timeStep',.0005,'speeds',[2.5,5,7.5,10],'dynamic',true,'initialSensorAmplitude',1e-5,'wakeLengthChords',1000,'coreRadius',1e-5,'modalCount',inf,'compactTangent',true,'staticRelaxation',.6,'checkStaticLimits',true);
names=fieldnames(settings);for k=1:numel(names),assert(isfield(o,names{k}),'cad:Options','Unknown FSI campaign setting.');o.(names{k})=settings.(names{k});end
saved=load(bundlePath,'bundle');b=saved.bundle;model=b.model;config=b.config;flow=b.flow;
if isfinite(o.modalCount),m=min(o.modalCount,size(model.K,1));model.K=model.K(1:m,1:m);model.M=model.M(1:m,1:m);model.C=model.C(1:m,1:m);model.PhiSurface=model.PhiSurface(:,1:m);model.observation=model.observation(:,1:m);end
config.discretization.dampingRatio=model.zeta;config.provenance.runtimeOverrides='assumed modal damping 0.008; literature proxy, not measured';
output=fullfile(root,'results','cad_fsi',char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS')));mkdir(output);cad.provenance(root,config,output);
staticRows=[];dynamicRows=[];dynamicResults=cell(numel(o.speeds),1);qPrevious=zeros(size(model.K,1),1);
for k=1:numel(o.speeds)
    flow.speed=o.speeds(k);folder=fullfile(output,"U_"+flow.speed);mkdir(folder);
    rigid=model.aero(model.X0,zeros(numel(model.X0)/3,3),flow);qOne=model.K\(model.PhiSurface'*reshape(rigid.nodalForces',[],1));
    equilibrium=fsi.static(model,flow,struct('initialQ',qPrevious,'relaxation',o.staticRelaxation,'loadTolerance',1e-7,'displacementTolerance',1e-7));qPrevious=equilibrium.q;
    yOne=model.observation*qOne;yTwo=model.observation*equilibrium.q;
    staticRows(end+1,:)=[flow.speed,rigid.CL,equilibrium.aero.CL,max(abs(yOne)),max(abs(yTwo)),equilibrium.relativeEquilibriumResidual,equilibrium.iterations]; %#ok<AGROW>
    writematrix(equilibrium.history,fullfile(folder,'static_residuals.csv'));save(fullfile(folder,'static.mat'),'equilibrium','qOne','rigid','-v7.3');
    if ~o.dynamic,continue;end
    % Consistent UVLM equilibrium, rather than imposing stationary VLM loads
    % as a time-periodic forcing. Tangent dynamics starts about this state.
    uvOptions=struct('wakeLength',o.wakeLengthChords*flow.referenceChord,'coreRadius',o.coreRadius);
    uvModel=model;uvModel.aero=@(X,V,f)fsi.uvlm_static(X,V,f,uvOptions);
    uvEquilibrium=fsi.static(uvModel,flow,struct('initialQ',equilibrium.q,'relaxation',o.staticRelaxation,'loadTolerance',1e-7,'displacementTolerance',1e-7));
    save(fullfile(folder,'uv_equilibrium.mat'),'uvEquilibrium','uvOptions','flow','-v7.3');
    X=reshape(model.X0,[],3)+reshape(model.PhiSurface*uvEquilibrium.q,3,[])';X=reshape(X,size(model.X0));
    nt=ceil(o.duration/o.timeStep)+1;time=linspace(0,o.duration,nt);dt=time(2)-time(1);
    fprintf('CAD FSI U=%g: UVLM tangent %dx%d panels, %d time steps\n',flow.speed,b.surface.nx,b.surface.ny,nt-1);
    fluid=uvlm.linear_model(X,model.PhiSurface,flow,dt,nt,struct('wakeLength',o.wakeLengthChords*flow.referenceChord,'coreRadius',o.coreRadius,'compactTangent',o.compactTangent));
    q0=zeros(size(model.K,1),1);q0(1)=o.initialSensorAmplitude/max(abs(model.observation(:,1)));v0=zeros(size(q0));
    save(fullfile(folder,'fluid_tangent.mat'),'fluid','model','-v7.3');
    dynamics=fsi.linear_dynamic(model,fluid,time,q0,v0);Y=model.observation*dynamics.Q;
    envelope=cad.decay_fit(time,dynamics.Q(1,:));
    envelope.observable='first structural modal coordinate; TE histories retained separately';metrics=ael.response_metrics(time,Y);
    assert(dynamics.relativeEnergyBalance<.02 && max(dynamics.relativeResidual)<1e-6,'cad:FSI','Dynamic energy/equilibrium check failed.');
    dynamicRows(end+1,:)=[flow.speed,envelope.frequencyHz,envelope.growthPerSecond,envelope.peakCount, ...
        max(metrics.fluctuationRMS),dynamics.relativeEnergyBalance,max(dynamics.relativeResidual)]; %#ok<AGROW>
    scope=struct('initialPerturbation_m',o.initialSensorAmplitude,'dampingRatioAssumed',model.zeta, ...
        'circulationRetained',true,'wakeModel','linearized UVLM, frozen ambient convection', ...
        'nonlinearUVLMComparison','independently verified for small perturbations; finite wake checked separately', ...
        'growthFitQualified',envelope.qualified,'experimentalValidation',false, ...
        'scope','conditional attached-flow small-perturbation FSI; no separation or viscous-deficit claim');
    cad.write_json(fullfile(folder,'dynamic_scope.json'),scope);cad.write_json(fullfile(folder,'growth_fit.json'),envelope);
    writematrix([time',Y'],fullfile(folder,'te_history.csv'));save(fullfile(folder,'dynamic.mat'),'dynamics','Y','uvEquilibrium','flow','scope','metrics','-v7.3');
    fig=figure('Visible','off','Color','w','Position',[100,100,1100,700]);if isprop(fig,'Theme'),fig.Theme='light';end
    tiledlayout(fig,2,1,'Padding','compact');ax=nexttile;plot(ax,time,1e3*Y);grid(ax,'on');xlabel(ax,'t [s]');ylabel(ax,'perturbación w_TE [mm]');title(ax,sprintf('CAD U=%g m/s, UVLM acoplado; amortiguamiento supuesto %.2f%%',flow.speed,100*model.zeta));
    ax=nexttile;plot(ax,time,dynamics.energy-dynamics.energy(1),time,dynamics.aerodynamicWork-dynamics.dissipation);grid(ax,'on');xlabel(ax,'t [s]');ylabel(ax,'energía [J]');legend(ax,'Cambio estructural','Trabajo aerodinámico - disipación','Location','best');
    exportgraphics(fig,fullfile(folder,'coupled_response.png'),'Resolution',150);close(fig);
    dynamicResults{k}=struct('time',time,'Y',Y,'growth',envelope,'relativeEnergyBalance',dynamics.relativeEnergyBalance);
    fprintf('CAD FSI U=%g: f=%.4g Hz, growth=%.4g /s, energy defect %.3g\n',flow.speed,envelope.frequencyHz,envelope.growthPerSecond,dynamics.relativeEnergyBalance);
end
static=array2table(staticRows,'VariableNames',{'U_ms','rigid_CL','coupled_CL','one_way_abs_TE_m','two_way_abs_TE_m','residual','iterations'});writetable(static,fullfile(output,'static.csv'));
dynamic=table();if ~isempty(dynamicRows),dynamic=array2table(dynamicRows,'VariableNames',{'U_ms','frequency_Hz','growth_per_s','peaks','maximum_RMS_m','relative_energy_defect','relative_equilibrium_residual'});writetable(dynamic,fullfile(output,'dynamic.csv'));end
% Limiting checks and tolerance/relaxation sensitivity use identical physics.
if o.checkStaticLimits
flow.speed=o.speeds(end);rigid=fsi.vlm(model.X0,zeros(numel(model.X0)/3,3),flow);stiff=model;stiff.K=model.K*1e6;
high=fsi.static(stiff,flow,struct('relaxation',.7,'loadTolerance',1e-7,'displacementTolerance',1e-7));
slow=fsi.static(model,flow,struct('relaxation',.35,'loadTolerance',1e-8,'displacementTolerance',1e-8));
normal=fsi.static(model,flow,struct('relaxation',.7,'loadTolerance',1e-7,'displacementTolerance',1e-7));
limits=struct('highStiffnessRelativeCL',abs(high.aero.CL/rigid.CL-1),'relaxationDisplacementRelative', ...
    norm(model.observation*(slow.q-normal.q),inf)/max(norm(model.observation*normal.q,inf),1e-15));
assert(limits.highStiffnessRelativeCL<1e-5 && limits.relaxationDisplacementRelative<1e-4,'cad:FSI','Static limit checks failed.');
else
 limits=struct('executedInThisRun',false,'scope','refinement replay; independent static limits and guide-flow relaxation checks must be cited from their own campaigns');
end
cad.write_json(fullfile(output,'static_limits.json'),limits);
report=struct('output',output,'static',static,'dynamic',dynamic,'limits',limits,'dynamicResults',{dynamicResults}, ...
    'scope','conditional CAD numerical campaign; refinement/uncertainty campaign required before final claims');
save(fullfile(output,'campaign.mat'),'report','o','-v7.3');fprintf('CAD FSI output: %s\n',output);
end
