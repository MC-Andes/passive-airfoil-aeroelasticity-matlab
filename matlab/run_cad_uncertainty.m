function report = run_cad_uncertainty(tangentPath)
%RUN_CAD_UNCERTAINTY Deterministic conditional ranges, never statistical CI.
% Uniform isotropic E and rho scale full-FEM K/M exactly. The aerodynamic
% tangent is reused around its declared baseline equilibrium; each static
% case is separately converged with the matching steady UVLM ring/wake model.
root=fileparts(fileparts(mfilename('fullpath')));s=load(tangentPath,'fluid','model');fluid=s.fluid;base=s.model;flow=fluid.flow;
uvOptions=fluid.settings;base.aero=@(X,V,f)fsi.uvlm_static(X,V,f,uvOptions);
parameters=[1,1,.008;.78,1,.008;3.05/2.617,1,.008;1,1124.6/1160,.008;1,1240/1160,.008;1,1,.004;1,1,.016;.78,1240/1160,.004;3.05/2.617,1124.6/1160,.016];rows=[];histories=cell(size(parameters,1),1);time=(0:fluid.wakeRows-2)*fluid.dt;
saved=load(fullfile(fileparts(tangentPath),'dynamic.mat'),'dynamics','uvEquilibrium');assert(numel(saved.dynamics.time)==numel(time) && max(abs(saved.dynamics.time-time))<1e-10,'cad:Time','Baseline response does not match this tangent.');
output=fullfile(fileparts(fileparts(tangentPath)),'uncertainty');if ~isfolder(output),mkdir(output);end
for k=1:size(parameters,1)
 model=base;model.K=base.K*parameters(k,1);model.M=base.M*parameters(k,2);omega=sqrt(diag(model.K)./diag(model.M));model.C=diag(2*parameters(k,3)*omega.*diag(model.M));
 initialQ=model.K\fluid.load0;
 q0=zeros(size(model.K,1),1);q0(1)=1e-5/max(abs(model.observation(:,1)));
 if k==1 && norm(model.C-base.C,'fro')<1e-12*norm(base.C,'fro') && norm(saved.dynamics.Q(:,1)-q0)<1e-14 && norm(saved.dynamics.V(:,1))<1e-14,equilibrium=saved.uvEquilibrium;dynamics=saved.dynamics;
 else,equilibrium=fsi.static(model,flow,struct('initialQ',initialQ,'relaxation',1,'loadTolerance',1e-7,'displacementTolerance',1e-7));dynamics=fsi.linear_dynamic(model,fluid,time,q0,zeros(size(q0)));end
 fit=cad.decay_fit(time,dynamics.Q(1,:));Y=model.observation*dynamics.Q;metrics=ael.response_metrics(time,Y);
 rows(end+1,:)=[k,2.617e9*parameters(k,1),1160*parameters(k,2),parameters(k,3),max(abs(model.observation*equilibrium.q)),equilibrium.aero.CL,omega(1)/(2*pi),fit.frequencyHz,fit.growthPerSecond,fit.peakCount,max(metrics.fluctuationRMS),dynamics.relativeEnergyBalance]; %#ok<AGROW>
 histories{k}=struct('time',time,'Y',Y,'growth',fit);fprintf('CAD uncertainty %d: E %.3g GPa rho %.4g zeta %.4g, growth %.4g\n',k,2.617*parameters(k,1),1160*parameters(k,2),parameters(k,3),fit.growthPerSecond);
end
report=array2table(rows,'VariableNames',{'case','E_Pa','rho_kg_m3','zeta_assumed','maximum_static_TE_m','CL','dry_frequency_Hz','coupled_frequency_Hz','growth_per_s','peaks','maximum_RMS_m','energy_defect'});writetable(report,fullfile(output,'material_damping.csv'));save(fullfile(output,'histories.mat'),'histories','parameters','-v7.3');
cad.write_json(fullfile(output,'scope.json'),struct('probabilityModel','none; engineering envelope, not confidence interval','staticModel','steady UVLM using the same wake/core as the dynamic baseline','sources',{{'Malim2023 E=2617MPa rho=1160; Morettini2022 damping=0.008 and density=1124.6; upper E 3.05GPa is an engineering sensitivity limit'}},'assumptions','22% lower effective E as orientation proxy; upper E is assumed; damping 0.004–0.016 is factor-two engineering sensitivity, not measured distribution','limitations','dynamic tangent fixed at baseline deformed equilibrium; not a new specimen calibration or a band for aerodynamic model discrepancy','sourceTangent',tangentPath));
fig=figure('Visible','off','Color','w','Position',[100,100,1100,700]);if isprop(fig,'Theme'),fig.Theme='light';end;tiledlayout(2,1);nexttile;plot(report.E_Pa/1e9,1e3*report.maximum_static_TE_m,'o');xlabel('E [GPa]');ylabel('máximo |w_TE| [mm]');title('Equilibrio UVLM, U=10 m/s');grid on;nexttile;scatter(report.zeta_assumed,report.growth_per_s,60,report.rho_kg_m3,'filled');xlabel('amortiguamiento supuesto');ylabel('crecimiento [1/s]');cb=colorbar;cb.Label.String='densidad [kg/m^3]';grid on;exportgraphics(fig,fullfile(output,'conditional_bands.png'),'Resolution',150);close(fig);
end
