function report = run_cad_viscous_postprocess(folders,outputFolder)
%RUN_CAD_VISCOUS_POSTPROCESS Common physical wake profiles and skin-offset audit.
root=fileparts(fileparts(mfilename('fullpath')));if nargin<2,outputFolder=fullfile(root,'results/cad_original/viscous_review');end
tab=readmatrix(fullfile(root,'data/reference/section_03.csv'));c=cad.load_case(root,fullfile(root,'cases/cad_original/extruded-fine.json'));profile=[(tab(:,2)-c.geometry.origin_m(1))/c.geometry.chord_m,(tab(:,3:4)-tab(1,5))/c.geometry.chord_m];span=c.geometry.span_m/c.geometry.chord_m;stations=[1.25,1.5,2,2.5];height=linspace(-.75,.75,241)';profiles=cell(numel(folders),1);rows=[];previous=[];
for k=1:numel(folders)
 s=load(fullfile(folders{k},'flow.mat'),'result','summary');r=s.result;r.summary=s.summary;g=r.grid;d=g.dimension;axes=cell(1,d);for j=1:d,index=repmat({1},1,d);index{j}=':';axes{j}=reshape(g.velocityCoordinates{1,j}(index{:}),[],1);end
 interpolation=griddedInterpolant(axes,r.meanVelocity{1},'linear','none');Y=zeros(numel(height),numel(stations));for j=1:numel(stations),if d==2,Y(:,j)=1-interpolation(stations(j)+zeros(size(height)),height);else,Y(:,j)=1-interpolation(stations(j)+zeros(size(height)),height,zeros(size(height)));end,end
 profiles{k}=struct('source',folders{k},'count',g.count,'stations',stations,'height',height,'deficit',Y,'Re',s.summary.Re,'dimension',d);change=NaN;if ~isempty(previous) && previous.dimension==d && previous.Re==s.summary.Re,change=norm(Y-previous.deficit,'fro')/max(norm(Y,'fro'),realmin);end;previous=profiles{k};rows(end+1,:)=[k,d,s.summary.Re,g.count(1),max(Y,[],'all'),change]; %#ok<AGROW>
 if d==2
  o=jsondecode(fileread(fullfile(folders{k},'case.json')));offsetRows=[];
  for offset=[1,1.5,2,3]
   skin=nsflow.surface_traction(r,profile,span,o.alpha,offset);offsetRows(end+1,:)=[offset,skin.forceCoefficient,s.summary.meanForceCoefficient,skin.pressureReferenceValue]; %#ok<AGROW>
   if offset==2
    writematrix(skin.rows,fullfile(folders{k},'surface_traction_upstream_reference.csv'));cad.write_json(fullfile(folders{k},'surface_reference.json'),rmfield(skin,'rows'));
    fig=figure('Visible','off','Color','w');if isprop(fig,'Theme'),fig.Theme='light';end;hold on;
    for side=1:2,take=skin.rows(:,1)==side;plot(skin.rows(take,2),skin.rows(take,7),'DisplayName',sprintf('cara %d',side));end
    set(gca,'YDir','reverse');grid on;xlabel('x/c');ylabel('Cp, referencia reservorio aguas arriba');legend;title(sprintf('CAD Re=%g, offset=%g c',s.summary.Re,skin.offsetChord));exportgraphics(fig,fullfile(folders{k},'surface_cp_upstream_reference.png'),'Resolution',150);close(fig);
   end
  end
  writematrix(offsetRows,fullfile(folders{k},'skin_offset_audit.csv'));
 end
end
report=array2table(rows,'VariableNames',{'case','dimension','Re','nx','maximum_deficit_on_profiles','relative_profile_change'});out=outputFolder;if ~isfolder(out),mkdir(out);end;writetable(report,fullfile(out,'wake_convergence.csv'));save(fullfile(out,'profiles.mat'),'profiles','report');
fig=figure('Visible','off','Color','w','Position',[100,100,1300,700]);if isprop(fig,'Theme'),fig.Theme='light';end;tiledlayout(1,numel(stations));for j=1:numel(stations),nexttile;hold on;for k=1:numel(profiles),plot(profiles{k}.deficit(:,j),height,'DisplayName',mat2str(profiles{k}.count));end;grid on;xlabel('1-u/U');ylabel('altura/c');title(sprintf('x/c=%g',stations(j)));end;legend;exportgraphics(fig,fullfile(out,'wake_profiles.png'),'Resolution',150);close(fig);disp(report);
end
