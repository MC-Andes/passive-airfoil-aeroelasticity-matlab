function report = run_cad_cavity_sensitivity(tangentPath)
%RUN_CAD_CAVITY_SENSITIVITY Conditional open-port impedance, fixed equilibrium.
% Static pressure is ambient in all three cases. The positive mass/damping
% proxies come from actual cavity volume gradients, not numerical springs.
root=fileparts(fileparts(mfilename('fullpath')));s=load(tangentPath,'fluid','model');p=load(fullfile(root,'results/cad_original/cavity_screen/operators.mat'),'M','C','scope');time=(0:s.fluid.wakeRows-2)*s.fluid.dt;rows=[];
saved=load(fullfile(fileparts(tangentPath),'dynamic.mat'),'dynamics');assert(numel(saved.dynamics.time)==numel(time) && max(abs(saved.dynamics.time-time))<1e-10,'cad:Time','Baseline response does not match this tangent.');
for factor=[0,1/3,1]
 m=s.model;m.M=m.M+factor*p.M;m.C=m.C+factor*p.C;q0=zeros(size(m.K,1),1);q0(1)=1e-5/max(abs(m.observation(:,1)));if factor==0 && norm(saved.dynamics.Q(:,1)-q0)<1e-14 && norm(saved.dynamics.V(:,1))<1e-14,r=saved.dynamics;else,r=fsi.linear_dynamic(m,s.fluid,time,q0,zeros(size(q0)));end;fit=cad.decay_fit(time,r.Q(1,:));rows(end+1,:)=[factor,fit.frequencyHz,fit.growthPerSecond,r.relativeEnergyBalance]; %#ok<AGROW>
end
report=array2table(rows,'VariableNames',{'vent_proxy_fraction','frequency_Hz','growth_per_s','energy_defect'});out=fullfile(root,'results/cad_original/cavity_screen');writetable(report,fullfile(out,'dynamic_sensitivity.csv'));cad.write_json(fullfile(out,'sensitivity_scope.json'),struct('ambientStaticPressure',true,'models','0: ambient dynamic pressure; 1/3: uniform-span source proxy; 1: lumped midspan source proxy','limitations','frequency-local oscillatory impedance estimate for open original ports; no closed-port prediction or measured flow resistance','source',tangentPath));disp(report);
end
