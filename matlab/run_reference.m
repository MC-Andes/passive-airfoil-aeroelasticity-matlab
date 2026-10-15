% RUN_REFERENCE CAD-derived surface with an EXPLICIT illustrative structure.
% It is NOT a prediction for the supplied manufactured airfoil.
projectRoot=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot,'matlab'));
outputFolder=fullfile(projectRoot,'results','reference');
if ~isfolder(outputFolder), mkdir(outputFolder); end
[surface,cadInfo]=pae.surface_from_sections(fullfile(projectRoot,'data','reference'),4,14,"constant");
% Replace this line with your identified/imported reduced structural model.
structure=pae.example_structure(surface);
model=pae.make_model(surface,structure);
flow=struct('rho',1.225,'speed',10,'alpha',deg2rad(3), ...
    'referencePoint',[surface.referenceChord/4,0,0]);
[summary,solutions]=pae.sweep(model,[0,5,10,15],deg2rad([0,3]),flow);
fig=pae.plot_study(model,summary,solutions);
writetable(summary,fullfile(outputFolder,'illustrative_sweep.csv'));
exportgraphics(fig,fullfile(outputFolder,'illustrative_study.png'),'Resolution',160);
save(fullfile(outputFolder,'illustrative_study.mat'),'model','summary','solutions','cadInfo');
disp(summary);
fprintf('CAD: span=%.6f m, mean chord=%.6f m, reference area=%.6f m^2.\n', ...
    surface.span,surface.referenceChord,surface.referenceArea);
fprintf('Estructura: %s\n',structure.source);
