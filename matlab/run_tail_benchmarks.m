function report = run_tail_benchmarks(refinements)
%RUN_TAIL_BENCHMARKS Reproducible comparison of nine published PLA plate cases.
% No parameter fitting; normalized frequencies and first-mode mesh changes.
if nargin<1, refinements=[8,16]; end
assert(numel(refinements)>=2 && all(diff(refinements)>0), ...
    'ael:Mesh','Supply at least two increasing mesh refinements.');
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'matlab'));
out=fullfile(root,'results','tail_benchmarks'); if ~isfolder(out), mkdir(out); end
cases=ael.reference_cases(); rows=cell(0,12);
for k=1:numel(cases)
    c=cases(k); previous=NaN;
    for level=1:numel(refinements)
        [s,fixed]=ael.build_plate_case(c,refinements(level));
        O=plate.probe(s.mesh,[c.chord,c.span],1);
        model=ael.reduce(s.K,s.M,fixed,6,.01,s.P,O);
        D=c.material.E*c.rigidThickness^3/(12*(1-c.material.nu^2));
        Omega=model.omega*c.span^2*sqrt(c.material.rho*c.rigidThickness/D);
        change=abs(Omega(1)-previous)/Omega(1); previous=Omega(1);
    end
    reference=c.reference;
    for rank=1:3
        target=reference.omega_nd_mean(rank); sd=reference.omega_nd_sample_sd(rank);
        rows(end+1,:)={c.id,rank,reference.source_mode_label(rank),Omega(rank), ...
            target,sd,100*(Omega(rank)/target-1),model.frequencyHz(rank), ...
            size(s.mesh.nodes,1),100*change,rank==1,string(c.sourceDOI)}; %#ok<AGROW>
    end
    fprintf('%s: Omega1 %.5f vs %.5f, error %+.2f%%, mesh change %.3f%%\n', ...
        c.id,Omega(1),reference.omega_nd_mean(1),100*(Omega(1)/reference.omega_nd_mean(1)-1),100*change);
end
report=cell2table(rows,'VariableNames',{'case_id','frequency_rank','source_mode_label', ...
    'model_Omega','experiment_Omega','experiment_sample_sd','error_percent', ...
    'nominal_model_Hz','nodes','first_mode_mesh_change_percent','first_bending_identified','source_DOI'});
writetable(report,fullfile(out,'modal_comparison.csv'));
first=report(report.frequency_rank==1,:);
fig=figure('Color','w','Visible','off','WindowStyle','normal');
if isprop(fig,'Theme'), fig.Theme='light'; end
errorbar(1:height(first),first.experiment_Omega,first.experiment_sample_sd,'o'); hold on;
plot(1:height(first),first.model_Omega,'x'); grid on;
xticks(1:height(first)); xticklabels(replace(erase(first.case_id,"HZ_"),"_"," "));
set(gca,'TickLabelInterpreter','none'); xtickangle(35); xlim([.5,height(first)+.5]);
ylabel('\Omega_1'); legend('Experimento ± SD','MITC4','Location','best');
title('Primer modo de flexión');
exportgraphics(fig,fullfile(out,'modal_comparison.png'),'Resolution',160); close(fig);
end
