function review = run_cad_final_review(manifestPath)
%RUN_CAD_FINAL_REVIEW Review explicit completed campaigns; never selects latest.
% Manifest paths: structural, geometryComparison, fullReducedComparison,
% loadComparison, partitionedComparison, panel96, panel128, timeHalf128,
% core, wake, speeds, viscous2D (three folders), viscous3D (two folders).
root=fileparts(fileparts(mfilename('fullpath')));p=jsondecode(fileread(manifestPath));
out=fullfile(root,'results/cad_original/final_review');if ~isfolder(out),mkdir(out);end
config=cad.load_case(root,fullfile(root,'cases/cad_original/extruded-fine.json'));cad.provenance(root,config,out);
structural=jsondecode(fileread(fullfile(p.structural,'acceptance.json')));
review=struct('createdUTC',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')), ...
 'manifest',p,'sourceStepSha256',config.geometry.sourceSha256,'experimentalValidation',false, ...
 'scope','conditional numerical development of the original connected CAD; assumed material/mount/guide; no specimen or tunnel validation');
review.structural=struct('numericalClosure',structural.numericalClosure,'discretizationCriteriaMet',structural.discretizationCriteriaMet, ...
 'algebraicCriteriaMet',structural.algebraicCriteriaMet,'relevantModeMask',structural.relevantModeMask, ...
 'historicalFlagExplanation','original acceptance.json structuralModelClosed=false and blockedStages are preserved; subsequent mechanism/flow evidence is reviewed here');
review.curvedGeometry=jsondecode(fileread(p.geometryComparison));
if isfield(p,'curvedQuadrature'),review.curvedQuadrature=jsondecode(fileread(p.curvedQuadrature));end
review.fullReduced=jsondecode(fileread(p.fullReducedComparison));review.fullLoad=jsondecode(fileread(p.loadComparison));
review.partitioned=jsondecode(fileread(p.partitionedComparison));
review.independentCriteria=struct('geometry',review.curvedGeometry.tailSensorRelative<.02 && max(review.curvedGeometry.frequencyRelative)<.01 && review.curvedGeometry.energyResidual<1e-6, ...
 'fullReduction',review.fullReduced.maximumPointwiseRelativeError<.02 && review.fullReduced.fullRelativeEnergyDefect<.02, ...
 'refinedLoads',max(review.fullLoad.fullReducedSensorRelative)<.02 && max(review.fullLoad.energyResidual)<1e-6, ...
 'partitioned',review.partitioned.nonlinearTangentSensorRelative<.02 && review.partitioned.relaxationToleranceSensorRelative<.02 && review.partitioned.energyDefect<.02 && review.partitioned.interfaceWorkDefect_J<1e-10);
if isfield(p,'curvedQuadrature'),review.independentCriteria.curvedQuadrature=review.curvedQuadrature.passed;end
if isfield(p,'pressureHistory')
 review.pressureHistory=jsondecode(fileread(p.pressureHistory));review.pressureObserver=jsondecode(fileread(p.pressureObserverComparison));review.independentCriteria.pressureHistory=review.pressureHistory.energyDefect<.02 && review.pressureHistory.interfaceWorkDefect_J<1e-10 && review.pressureObserver.passed;
end
if isfield(p,'posteriorScope')
 review.posterior=jsondecode(fileread(p.posteriorScope));t=readtable(p.posteriorConvergence);review.posterior.curvatureFieldLastRelativeChange=t.relative_curvature_field_change(end);review.independentCriteria.posteriorAnalytical=review.posterior.analyticalParabolaMaximumError_per_m<1e-9;
end
if isfield(p,'frequencyBenchmarks')
 benchmarks=cell(numel(p.frequencyBenchmarks),1);for k=1:numel(benchmarks),t=readtable(p.frequencyBenchmarks{k});t=t(end-1:end,:);benchmarks{k}=struct('source',p.frequencyBenchmarks{k},'refinedCases',table2struct(t),'passed',all(t.amplitude_relative_error<.02 & t.phase_error_deg<2));end
 review.uvlmBenchmarks=benchmarks;review.independentCriteria.unsteadyAnalytical=all(cellfun(@(x)x.passed,benchmarks));
end
review.static=struct();staticNames={'staticHorseshoe','staticUVLM'};
for k=1:numel(staticNames)
 name=staticNames{k};if ~isfield(p,name),continue;end;t=readtable(p.(name));changes=t.sensor_vector_relative_change;take=changes(max(2,end-1):end);review.static.(name)=struct('cases',table2struct(t),'lastTwoChanges',take,'passed',height(t)>=3 && all(take<.02));
end
names={'panel96','panel128','timeHalf128','core','wake'};rows=[];
for k=1:numel(names)
 folder=p.(names{k});g=jsondecode(fileread(fullfile(folder,'U_10/growth_fit.json')));t=readtable(fullfile(folder,'dynamic.csv'));
 assert(g.qualified,'cad:Review','Growth fit has too few qualified peaks: %s',folder);
 rows(end+1,:)=[k,g.frequencyHz,g.growthPerSecond,g.peakCount,max(t.relative_energy_defect),max(t.relative_equilibrium_residual)]; %#ok<AGROW>
end
dynamic=array2table(rows,'VariableNames',{'case','frequency_Hz','growth_per_s','peaks','energy_defect','equilibrium_residual'});dynamic.case_name=string(names(:));writetable(dynamic,fullfile(out,'dynamic_cases.csv'));
timeCoarse=p.panel128;if isfield(p,'timeCoarse128'),timeCoarse=p.timeCoarse128;end
space=cad.compare_free_response(p.panel96,p.panel128);time=cad.compare_free_response(timeCoarse,p.timeHalf128);core=cad.compare_free_response(p.panel128,p.core);wake=cad.compare_free_response(p.panel128,p.wake);
review.dynamic=struct('space',space,'time',time,'core',core,'wake',wake,'allGrowthFitsQualified',true, ...
 'numericalCriteriaMet',space.passed && time.passed && core.passed && wake.passed, ...
 'criteria','frequency 1%, growth/RMS/maximum history amplitude normalized by reference peak 2%, phase 2 degrees; fitted decaying-envelope amplitude is diagnostic; pairwise changes are not certified total error bounds', ...
 'flutterEstablished',false);
if isfield(p,'speedCases')
 speeds=table();for k=1:numel(p.speedCases),speeds=[speeds;cad.response_row(p.speedCases(k).folder,p.speedCases(k).speed)];end %#ok<AGROW>
else
 speeds=readtable(fullfile(p.speeds,'dynamic.csv'));
end
reference=cad.response_row(p.timeHalf128,10);speeds=[speeds;reference];writetable(speeds,fullfile(out,'speed_map.csv'));review.speedMap=table2struct(speeds);
if isfield(p,'speeds96')
 speedSpace=cell(3,1);for k=1:3,if isfield(p,'speedSpatialCases'),folder=p.speedSpatialCases(k).folder;elseif isfield(p,'speedCases'),assert(p.speedCases(k).speed==2.5*k,'cad:Review','Speed cases must be ordered 2.5, 5, 7.5 m/s.');folder=p.speedCases(k).folder;else,folder=p.speeds;end;speedSpace{k}=cad.compare_free_response(p.speeds96,folder,2.5*k);end
 review.dynamic.speedSpace=speedSpace;review.dynamic.speedSpatialCriteriaMet=all(cellfun(@(x)x.passed,speedSpace));
 review.dynamic.numericalCriteriaMet=review.dynamic.numericalCriteriaMet && review.dynamic.speedSpatialCriteriaMet;
end
if isfield(p,'speedTimePairs')
 speedTime=cell(numel(p.speedTimePairs),1);for k=1:numel(speedTime),pair=p.speedTimePairs(k);speedTime{k}=cad.compare_free_response(pair.coarse,pair.fine,pair.speed);end
 review.dynamic.speedTime=speedTime;review.dynamic.speedTemporalCriteriaMet=all(cellfun(@(x)x.passed,speedTime));review.dynamic.numericalCriteriaMet=review.dynamic.numericalCriteriaMet && review.dynamic.speedTemporalCriteriaMet;
end
if isfield(p,'guides')
 guideNames={'B','C'};for k=1:2,name=guideNames{k};comparison=cad.compare_free_response(p.guides.([name,'_coarse']),p.guides.([name,'_fine']));review.guideDynamic.(name)=comparison;review.dynamic.numericalCriteriaMet=review.dynamic.numericalCriteriaMet && comparison.passed;end
end
vrows=[];for k=1:numel(p.viscous2D),s=jsondecode(fileread(fullfile(p.viscous2D{k},'summary.json')));vrows(end+1,:)=[s.count(:)',s.meanForceCoefficient(:)',s.maximumDivergence,s.maximumMomentumDefect];end %#ok<AGROW>
vt=array2table(vrows,'VariableNames',{'nx','ny','CD','CL','maximum_divergence','momentum_defect'});writetable(vt,fullfile(out,'viscous_2d.csv'));
review.viscous=struct('lastRelativeCDChange',abs(vt.CD(end)/vt.CD(end-1)-1),'lastRelativeCLChange',abs(vt.CL(end)/vt.CL(end-1)-1), ...
 'scope','laminar CAD exterior, Re100 2D and Re20 capped 3D; pressure/force/wake accuracy assessed separately; no tunnel-Re extrapolation');
three=cell(numel(p.viscous3D),1);for k=1:numel(three),three{k}=jsondecode(fileread(fullfile(p.viscous3D{k},'summary.json')));end;review.viscous.threeDimensional=three;
review.verification=jsondecode(fileread(fullfile(root,'results/verification/summary.json')));
review.acceptedConditionalNumericalModel=review.structural.numericalClosure && review.dynamic.numericalCriteriaMet && review.verification.passed;
fields=fieldnames(review.independentCriteria);for k=1:numel(fields),review.acceptedConditionalNumericalModel=review.acceptedConditionalNumericalModel && review.independentCriteria.(fields{k});end
for k=1:numel(staticNames),name=staticNames{k};if isfield(review.static,name),review.acceptedConditionalNumericalModel=review.acceptedConditionalNumericalModel && review.static.(name).passed;end,end
review.notEstablished={'specimen material/print/mount identification','actual slider clearance/friction/preload','tunnel Cp and turbulent wake','experimental flutter boundary','PIV/anemometry comparison'};
cad.write_json(fullfile(out,'acceptance.json'),review);copyfile(manifestPath,fullfile(out,'campaign_inputs.json'));
files={manifestPath,fullfile(p.structural,'acceptance.json'),p.geometryComparison,p.fullReducedComparison,p.loadComparison,p.partitionedComparison};
for k=1:numel(names),files{end+1}=fullfile(p.(names{k}),'dynamic.csv');files{end+1}=fullfile(p.(names{k}),'U_10/growth_fit.json');end %#ok<AGROW>
if isfield(p,'timeCoarse128'),files{end+1}=fullfile(p.timeCoarse128,'dynamic.csv');files{end+1}=fullfile(p.timeCoarse128,'U_10/growth_fit.json');end
fields=fieldnames(p);for k=1:numel(fields),value=p.(fields{k});if ischar(value) && isfile(value),files{end+1}=value;end,end %#ok<AGROW>
if isfield(p,'frequencyBenchmarks'),files=[files,p.frequencyBenchmarks(:)'];end
folders=[p.viscous2D(:);p.viscous3D(:)];for k=1:numel(folders),files{end+1}=fullfile(folders{k},'summary.json');files{end+1}=fullfile(folders{k},'case.json');end %#ok<AGROW>
if isfield(p,'guides'),fields=fieldnames(p.guides);for k=1:numel(fields),files{end+1}=fullfile(p.guides.(fields{k}),'dynamic.csv');files{end+1}=fullfile(p.guides.(fields{k}),'U_10/growth_fit.json');end,end %#ok<AGROW>
if isfield(p,'speedCases'),for k=1:numel(p.speedCases),subfolder=fullfile(p.speedCases(k).folder,"U_"+p.speedCases(k).speed);files{end+1}=fullfile(subfolder,'growth_fit.json');files{end+1}=fullfile(subfolder,'dynamic.mat');end,end %#ok<AGROW>
if isfield(p,'speedSpatialCases'),for k=1:numel(p.speedSpatialCases),files{end+1}=fullfile(p.speedSpatialCases(k).folder,"U_"+p.speedSpatialCases(k).speed,'growth_fit.json');end,end %#ok<AGROW>
if isfield(p,'speedTimePairs'),for k=1:numel(p.speedTimePairs),files{end+1}=fullfile(p.speedTimePairs(k).coarse,"U_"+p.speedTimePairs(k).speed,'growth_fit.json');end,end %#ok<AGROW>
files=cellstr(unique(string(files),'stable'));
hashes=strings(numel(files),1);for k=1:numel(files),hashes(k)=cad.sha256(files{k});end;writetable(table(string(files(:)),hashes,'VariableNames',{'file','sha256'}),fullfile(out,'evidence_hashes.csv'));
disp(dynamic);disp(review.dynamic);
end
