function report = run_cad_modal_response(solutionFolder)
%RUN_CAD_MODAL_RESPONSE Conditional free decay of an observed CAD structural mode.
% No aerodynamic forcing. Modes come from a saved CAD campaign; their residuals
% and spatial convergence remain independent requirements. Damping is assumed.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'matlab'));
s=load(fullfile(solutionFolder,'solution.mat'),'config','modes','sensor','result');
config=s.config; modes=s.modes; O=s.sensor.O*modes.Phi; n=numel(modes.omega); z=config.discretization.dampingRatio;
assert(z>=0 && z<1 && max(abs(O(:,1)))>0,'cad:Dynamics','Require an observed underdamped mode.');
model=struct('K',diag(modes.omega.^2),'M',eye(n),'C',diag(2*z*modes.omega), ...
    'B',modes.Phi'*s.result.F,'O',O);
q0=zeros(n,1); q0(1)=1e-5/max(abs(O(:,1))); T=1/modes.frequencyHz(1); err=zeros(2,1); energyDefect=err;
for k=1:2
    time=linspace(0,10*T,1000*2^(k-1)+1); r=ael.transient(model,time,@(t)0,q0,zeros(n,1));
    wd=modes.omega(1)*sqrt(1-z^2);
    exact=q0(1)*exp(-z*modes.omega(1)*time).*(cos(wd*time)+z/sqrt(1-z^2)*sin(wd*time));
    err(k)=max(abs(r.Y-O(:,1)*exact),[],'all')/1e-5;
    dissipated=trapz(time,sum(r.V.*(model.C*r.V),1));
    energyDefect(k)=abs(r.energy(end)-r.energy(1)+dissipated)/r.energy(1);
end
assert(err(2)<.02 && err(2)<err(1) && energyDefect(2)<.02,'cad:Dynamics','Time convergence/energy check failed.');
out=fullfile(solutionFolder,'free_decay'); if ~isfolder(out), mkdir(out); end
writetable(array2table([time',r.Y',r.energy']),fullfile(out,'signals.csv'));
metrics=ael.response_metrics(time,r.Y); cad.write_json(fullfile(out,'metrics.json'),metrics);
report=struct('scope','prescribed initial perturbation; free structural modal decay, no flow excitation', ...
    'dampingStatus','assumed modal ratio, not measured','zeta',z,'initialMaximumSensor_m',1e-5, ...
    'timeWindow_s',[time(1),time(end)],'steps',[1000,2000],'maximumAmplitudeRelativeError',err, ...
    'relativeEnergyBalanceDefect',energyDefect,'maximumSourceEigenResidual',max(modes.eigenResidual), ...
    'spatialConvergenceValidated',false,'interpretation','temporal check of a conditional reduced model; not coupled CAD validation');
cad.write_json(fullfile(out,'scope.json'),report);
fig=figure('Visible','off','Color','w','Position',[100,100,1000,600]);
if isprop(fig,'Theme'), fig.Theme='light'; end
tiledlayout(fig,2,1,'Padding','compact'); ax=nexttile;
plot(ax,time,1e3*r.Y,'LineWidth',1); grid(ax,'on'); xlabel(ax,'t [s]'); ylabel(ax,'w_z cerca de TE [mm]');
title(ax,'Decaimiento modal libre: perturbación inicial impuesta, sin carga de flujo');
ax=nexttile; plot(ax,time,r.energy/r.energy(1),'LineWidth',1.2); grid(ax,'on');
xlabel(ax,'t [s]'); ylabel(ax,'E/E(0)'); title(ax,sprintf('Amortiguamiento modal supuesto: %.1f%%',100*z));
exportgraphics(fig,fullfile(out,'free_decay.png'),'Resolution',150); close(fig);
end
