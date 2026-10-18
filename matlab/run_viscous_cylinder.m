function report = run_viscous_cylinder(counts,tEnd)
%RUN_VISCOUS_CYLINDER Independent external-flow reference at Re=20.
% Dennis & Chang 1970 Cd=2.05 is a comparison target, not enforced/calibrated.
if nargin<1,counts=[192,256,384];end;if nargin<2,tEnd=30;end
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'matlab'));
output=fullfile(root,'results','viscous_cylinder',char(datetime('now','Format','yyyyMMdd_HHmmss')));mkdir(output);
rows=[];
for n=counts
    folder=fullfile(output,"n_"+n);h=16/n;mask=@(coords,t) cylinder_mask(coords,h/2);
    result=nsflow.external([n,n],[-6,10;-8,8],20,mask,tEnd,folder);
    Cd=result.summary.meanForceCoefficient(1);rows(end+1,:)=[n,Cd,abs(Cd/2.05-1),result.summary.maximumDivergence]; %#ok<AGROW>
    nsflow.plot_wake(result,fullfile(folder,'wake.png'),sprintf('Cilindro Re=20, %dx%d; dominio y penalización declarados',n,n));
end
report=array2table(rows,'VariableNames',{'n','mean_Cd','relative_error_to_DennisChang','maximum_divergence'});writetable(report,fullfile(output,'comparison.csv'));disp(report);
end
function mask=cylinder_mask(coords,width)
s=sqrt(coords{1}.^2+coords{2}.^2)-.5;mask=zeros(size(s));mask(s<=-width)=1;take=abs(s)<width;mask(take)=.5-.5*sin(pi*s(take)/(2*width));
end
