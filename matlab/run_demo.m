% RUN_DEMO Generalized passive aeroelastic example. Parameters are illustrative.
projectRoot=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot,'matlab'));
outputFolder=fullfile(projectRoot,'results','demo');
if ~isfolder(outputFolder), mkdir(outputFolder); end
[model,flow]=pae.demo_case(4,12);
[summary,solutions]=pae.sweep(model,[0,5,10,15,20],deg2rad([2,5]),flow);
equilibrium=pae.solve_static(model,flow);
localDiagnostic=pae.linearize(model,equilibrium);
% Slower aerodynamic variation than the dry structural modes. This example
% illustrates the integration interface; it is not a flutter simulation.
flowFunction=@(t) struct('rho',flow.rho,'speed',flow.speed, ...
    'alpha',flow.alpha+deg2rad(1)*sin(2*pi*.8*t),'referencePoint',flow.referencePoint);
response=pae.simulate(model,flowFunction,(0:.002:1)',equilibrium.q,[]);
fig=pae.plot_study(model,summary,solutions,response);
writetable(summary,fullfile(outputFolder,'static_sweep.csv'));
exportgraphics(fig,fullfile(outputFolder,'study.png'),'Resolution',160);
savefig(fig,fullfile(outputFolder,'study.fig'));
save(fullfile(outputFolder,'study.mat'),'model','flow','summary','solutions','response','localDiagnostic');
disp(summary); fprintf('Resultados del ejemplo: %s\n',outputFolder);
