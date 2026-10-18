function report = run_morphing_study(casePath)
%RUN_MORPHING_STUDY Generate explicit CFD plan; no unexecuted result values.
root=fileparts(fileparts(mfilename('fullpath')));if nargin<1,casePath=fullfile(root,'cases/morphing_wing/case.json');end
c=morphing.load_case(root,casePath);label=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS'));out=fullfile(root,'results/morphing_wing',['study_',label]);mkdir(out);
cad.provenance(root,c,out);rows=[];
for U=c.flow.speeds_m_s(:)',for a=c.flow.alpha_deg(:)',for j=0:c.slider.positions-1
 fraction=j/(c.slider.positions-1);r=morphing.aerodynamic_coefficients(zeros(0,3),zeros(0,3),c,U,a);
 rows(end+1,:)=[U,a,j,fraction,r.q_Pa,r.Re,1000*fraction*c.slider.referenceFinal_m(:)']; %#ok<AGROW>
end,end,end
t=array2table(rows,'VariableNames',{'U_ms','alpha_deg','position_index','fraction','q_Pa','Re','command_if_minus12p5_mm','command_if_minus25_mm'});
t.status=repmat("planned_pending_final_travel_and_transition_CFD",height(t),1);writetable(t,fullfile(out,'planned_cases.csv'));
report=struct('output',out,'caseCount',height(t),'geometrySha256',c.geometry.sourceSha256,'caseSha256',c.caseSha256,'primary',c.analysis.primary, ...
 'complementary',{c.analysis.complementary},'tunnelCFDExecuted',false,'finalSliderTravelConfirmed',false, ...
 'referenceEndpoints_m',c.slider.referenceFinal_m,'referenceArea_m2',c.geometry.referenceArea_m2,'referenceChord_m',c.geometry.chord_m,'airDensity_kg_m3',c.flow.rho_kg_m3, ...
 'readiness','geometry/configuration and numerical contracts updated; planned rows are not computed aerodynamic outcomes');
cad.write_json(fullfile(out,'summary.json'),report);disp(report);
end
