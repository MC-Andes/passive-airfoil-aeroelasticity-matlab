function report = run_cad_time_refinement(tangentPath)
%RUN_CAD_TIME_REFINEMENT Half dt, identical trim, exact independent wake rows.
root=fileparts(fileparts(mfilename('fullpath')));parent=load(tangentPath,'fluid','model');previousFolder=fileparts(fileparts(tangentPath));d=load(fullfile(fileparts(tangentPath),'dynamic.mat'),'dynamics','uvEquilibrium','flow');
provenance=jsondecode(fileread(fullfile(previousFolder,'provenance.json')));config=provenance.case;config.discretization.timeStep=parent.fluid.dt/2;
output=fullfile(root,'results/cad_fsi',char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS')));mkdir(output);cad.provenance(root,config,output);
cad.write_json(fullfile(output,'parent.json'),struct('tangent',tangentPath,'tangentSha256',cad.sha256(tangentPath),'dynamicSha256',cad.sha256(fullfile(fileparts(tangentPath),'dynamic.mat')),'method','exact vortex subdivision identities; steady trim reused because it is independent of dt'));
fluid=uvlm.refine_time_model(parent.fluid);model=parent.model;flow=d.flow;uvEquilibrium=d.uvEquilibrium;time=(0:fluid.wakeRows-2)*fluid.dt;
folder=fullfile(output,"U_"+flow.speed);mkdir(folder);save(fullfile(folder,'fluid_tangent.mat'),'fluid','model','-v7.3');
dynamics=fsi.linear_dynamic(model,fluid,time,d.dynamics.Q(:,1),d.dynamics.V(:,1));Y=model.observation*dynamics.Q;metrics=ael.response_metrics(time,Y);fit=cad.decay_fit(time,dynamics.Q(1,:));
assert(fit.qualified && dynamics.relativeEnergyBalance<.02 && max(dynamics.relativeResidual)<1e-6,'cad:Refine','Refined response failed its numerical checks.');
scope=struct('dampingRatioAssumed',model.zeta,'initialPerturbation_m',max(abs(Y(:,1))),'circulationRetained',true,'experimentalValidation',false,'wakeModel','frozen ambient convection, independent fine histories','scope','same CAD/trim/basis; exact half-step kernel construction independently verified against direct assembly');
save(fullfile(folder,'dynamic.mat'),'dynamics','Y','uvEquilibrium','flow','scope','metrics','-v7.3');cad.write_json(fullfile(folder,'growth_fit.json'),fit);cad.write_json(fullfile(folder,'dynamic_scope.json'),scope);writematrix([time',Y'],fullfile(folder,'te_history.csv'));
dynamic=cad.response_row(output,flow.speed);writetable(dynamic,fullfile(output,'dynamic.csv'));
if isfile(fullfile(previousFolder,'static.csv')),copyfile(fullfile(previousFolder,'static.csv'),fullfile(output,'static.csv'));
else
 old=load(fullfile(fileparts(tangentPath),'static.mat'),'equilibrium','qOne','rigid');static=array2table([flow.speed,old.rigid.CL,old.equilibrium.aero.CL,max(abs(model.observation*old.qOne)),max(abs(model.observation*old.equilibrium.q)),old.equilibrium.relativeEquilibriumResidual,old.equilibrium.iterations],'VariableNames',{'U_ms','rigid_CL','coupled_CL','one_way_abs_TE_m','two_way_abs_TE_m','residual','iterations'});writetable(static,fullfile(output,'static.csv'));
end
fig=figure('Visible','off','Color','w','Position',[100,100,1100,700]);if isprop(fig,'Theme'),fig.Theme='light';end;tiledlayout(2,1);nexttile;plot(time,1e3*Y);grid on;xlabel('t [s]');ylabel('perturbación w_TE [mm]');title(sprintf('CAD UVLM U=%g m/s, dt=%g s',flow.speed,fluid.dt));nexttile;plot(time,dynamics.energy-dynamics.energy(1),time,dynamics.aerodynamicWork-dynamics.dissipation);grid on;xlabel('t [s]');ylabel('energía [J]');legend('Cambio estructural','Trabajo aerodinámico - disipación');exportgraphics(fig,fullfile(folder,'coupled_response.png'),'Resolution',150);close(fig);
report=struct('output',output,'dynamic',dynamic,'parentTangent',tangentPath,'timeStep',fluid.dt,'scope',scope);save(fullfile(output,'campaign.mat'),'report','-v7.3');fprintf('CAD exact time refinement output: %s\n',output);disp(dynamic);
end
