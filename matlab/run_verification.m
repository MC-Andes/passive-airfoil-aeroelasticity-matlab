function report = run_verification()
%RUN_VERIFICATION Execute both physics verification suites and save reports.
% From the project root: addpath('matlab'); report=run_verification();
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'matlab'),fullfile(root,'matlab','tests'));
output=fullfile(root,'results','verification');
if ~isfolder(output), mkdir(output); end
aerodynamics=run_tests();
fem=run_fem_tests();
writetable(aerodynamics,fullfile(output,'aerodynamics_and_reduced_model.csv'));
writetable(fem,fullfile(output,'solid_fem.csv'));
report=struct('matlabVersion',version,'matlabRelease',version('-release'), ...
    'aerodynamicsChecks',height(aerodynamics),'femChecks',height(fem), ...
    'passed',all(aerodynamics.passed) && all(fem.passed), ...
    'scope','Numerical verification; not experimental validation or STEP mesh convergence');
fid=fopen(fullfile(output,'summary.json'),'w');
assert(fid~=-1,'pae:Report','Cannot create verification report.');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
assert(report.passed,'pae:Verification','Verification failed.');
fprintf('Verification passed: %d reduced/aerodynamic + %d solid FEM checks.\n', ...
    report.aerodynamicsChecks,report.femChecks);
end
