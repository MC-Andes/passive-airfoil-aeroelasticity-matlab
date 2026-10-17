function report = run_verification()
%RUN_VERIFICATION Execute all seven physics verification suites and save reports.
% From the project root: addpath('matlab'); report=run_verification();
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'matlab'),fullfile(root,'matlab','tests'));
output=fullfile(root,'results','verification');
if ~isfolder(output), mkdir(output); end
aerodynamics=run_tests();
fem=run_fem_tests();
tail=run_tail_tests();
cad=run_cad_tests();
mechanics=run_mechanics_extension_tests();
uvlm=run_uvlm_tests();
viscous=run_viscous_tests();
writetable(mechanics,fullfile(output,'mechanics_extensions.csv'));
writetable(uvlm,fullfile(output,'uvlm_fsi.csv'));
writetable(viscous,fullfile(output,'viscous.csv'));
writetable(aerodynamics,fullfile(output,'aerodynamics_and_reduced_model.csv'));
writetable(fem,fullfile(output,'solid_fem.csv'));
writetable(tail,fullfile(output,'trailing_edge.csv'));
writetable(cad,fullfile(output,'cad_original.csv'));
report=struct('matlabVersion',version,'matlabRelease',version('-release'), ...
    'aerodynamicsChecks',height(aerodynamics),'femChecks',height(fem), ...
    'trailingEdgeChecks',height(tail), ...
    'cadChecks',height(cad),'mechanicsChecks',height(mechanics),'uvlmChecks',height(uvlm),'viscousChecks',height(viscous), ...
    'passed',all(mechanics.passed) && all(uvlm.passed) && all(viscous.passed) && all(cad.passed) && all(aerodynamics.passed) && all(fem.passed) && all(tail.passed), ...
    'scope','Numerical verification; not experimental validation or STEP mesh convergence');
fid=fopen(fullfile(output,'summary.json'),'w');
assert(fid~=-1,'pae:Report','Cannot create verification report.');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
assert(report.passed,'pae:Verification','Verification failed.');
fprintf('Verification passed: %d total checks across seven physics suites.\n',height(aerodynamics)+height(fem)+height(tail)+height(cad)+height(mechanics)+height(uvlm)+height(viscous));
end
