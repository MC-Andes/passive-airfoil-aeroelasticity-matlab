function report = run_cad_viscous(settings)
%RUN_CAD_VISCOUS Own laminar NS on the actual CAD exterior section/extrusion.
% The Reynolds number is deliberately declared; this is not turbulent tunnel
% flow and no 2D result is used to claim finite-span tip-flow accuracy.
if nargin<1,settings=struct();end
root=fileparts(fileparts(mfilename('fullpath')));o=struct('count',[512,256],'bounds',[-2,5;-1.75,1.75],'Re',100,'duration',10,'alpha',deg2rad(3),'amplitude',0,'omega',2,'penaltyFactor',.05,'penaltyProjection','coupled','timeStepScale',1,'restartFolder','');
keys=fieldnames(settings);for k=1:numel(keys),assert(isfield(o,keys{k}));o.(keys{k})=settings.(keys{k});end
c=jsondecode(fileread(fullfile(root,'cases/cad_original/extruded-fine.json')));table=readmatrix(fullfile(root,'data/reference/section_03.csv'));chord=c.geometry.chord_m;
profile=[(table(:,2)-c.geometry.origin_m(1))/chord,(table(:,3:4)-table(1,5))/chord];
span=c.geometry.span_m/chord;o.smoothingWidth=min(diff(o.bounds,1,2)./o.count');
output=fullfile(root,'results','cad_viscous',char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS')));mkdir(output);
c.flow.speed_m_s=o.Re*c.flow.mu_Pa_s/(c.flow.rho_kg_m3*chord);c.flow.alpha_deg=rad2deg(o.alpha);c.flow.Re=o.Re;c.provenance.runtimeOverrides='laminar viscous branch; physical speed derived from the declared Reynolds number';cad.provenance(root,c,output);
mask=@(coords,t)nsflow.cad_mask(coords,t,profile,span,o);options=struct('penaltyFactor',o.penaltyFactor,'penaltyProjection',o.penaltyProjection,'timeStepScale',o.timeStepScale);
if numel(o.count)==3,options.referenceArea=span;end
if ~isempty(o.restartFolder),old=load(fullfile(o.restartFolder,'flow.mat'),'result');assert(isequal(old.result.grid.count,o.count) && isequal(old.result.grid.bounds,o.bounds),'nsflow:Restart','Restart grid mismatch.');options.initialState=old.result.state;restartSha=cad.sha256(fullfile(o.restartFolder,'flow.mat'));cad.write_json(fullfile(output,'restart_provenance.json'),struct('folder',o.restartFolder,'flowSha256',restartSha,'initialTime',old.result.state.time));end
if o.amplitude~=0,options.wallVelocity=@(coords,t,a)nsflow.cad_wall_velocity(coords,t,a,o);end
result=nsflow.external(o.count,o.bounds,o.Re,mask,o.duration,output,options);nsflow.plot_wake(result,output);
if numel(o.count)==3,nsflow.plot_cross_wake(result,output);end
if numel(o.count)==2
 surface=nsflow.surface_traction(result,profile,span,o.alpha,2);writematrix(surface.rows,fullfile(output,'surface_traction.csv'));cad.write_json(fullfile(output,'surface_balance.json'),rmfield(surface,'rows'));
 fig=figure('Visible','off','Color','w');if isprop(fig,'Theme'),fig.Theme='light';end;hold on;for side=1:2,take=surface.rows(:,1)==side;plot(surface.rows(take,2),surface.rows(take,7),'DisplayName',sprintf('cara %d',side));end;set(gca,'YDir','reverse');grid on;xlabel('x/c');ylabel('Cp muestreado fuera de capa penalizada');legend;title(sprintf('CAD laminar Re=%g, offset=%g c',o.Re,surface.offsetChord));exportgraphics(fig,fullfile(output,'surface_cp.png'),'Resolution',150);close(fig);
end
cad.write_json(fullfile(output,'case.json'),o);cad.write_json(fullfile(output,'geometry_scope.json'),struct('sourceSha256',c.geometry.sourceSha256,'profile','audited midspan CAD exterior envelope; data/reference/section_03.csv','spanOverChord',span,'internalCavities','sealed for exterior flow; not cavity-pressure prediction','spanEnds','capped in 3D; conditional hypothesis','incidenceConvention','positive flow incidence = negative body rotation in uniform horizontal flow','pressureGauge','domain pressure zero mean; skin Cp subtracts upstream reservoir pressure, recorded in surface_balance.json','wallMotion','prescribed heave only if amplitude is nonzero; not self-excited FSI','experimentalValidation',false));
report=struct('output',output,'summary',result.summary);fprintf('CAD viscous output: %s\n',output);
end
