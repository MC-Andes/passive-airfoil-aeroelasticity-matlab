function report = run_cad_posterior(campaign,dynamicFolder,basisFile)
%RUN_CAD_POSTERIOR Material posterior shapes, slopes and fitted curvatures.
root=fileparts(fileparts(mfilename('fullpath')));out=fullfile(root,'results/cad_original/posterior');if ~isfolder(out),mkdir(out);end
xi=.68:.005:.995;previous=[];rows=[];sample=[];
for level=1:3
 s=load(fullfile(campaign,"tet10_level_"+level,'solution.mat'),'config','result');c=s.config;if level==1,cad.provenance(root,c,out);end
 mesh=fem3d.elevate(fem3d.read_mesh(fullfile(root,c.discretization.meshes{level}),c.geometry.axisOrder,c.geometry.origin_m));eta=c.sensors.eta(:);H=sparse(0,size(mesh.nodes,1));points=[];
 for i=1:numel(xi),cc=c;cc.sensors.x_over_c=xi(i);sensor=cad.sensors(mesh,cc);H=[H;sensor.H];points=[points;sensor.points];end %#ok<AGROW>
 U=H*reshape(100*s.result.u,3,[])';metrics=cad.posterior_metrics(points,U,xi,eta);change=NaN;if ~isempty(previous),change=norm(metrics.curvatureChange_per_m-previous,'fro')/max(norm(metrics.curvatureChange_per_m,'fro'),realmin);end;previous=metrics.curvatureChange_per_m;
 rows(end+1,:)=[level,max(abs(U(:,3))),max(abs(metrics.slopeChange_rad),[],'all'),max(abs(metrics.curvatureChange_per_m),[],'all'),change]; %#ok<AGROW>
 writetable(array2table(metrics.rows,'VariableNames',{'x_over_c','eta','slope_change_rad','curvature_change_per_m','reference_curvature_per_m','deformed_curvature_per_m'}),fullfile(out,"pressure_curvature_level_"+level+".csv"));
 if level==3,sample=struct('H',H,'points',points,'xi',xi,'eta',eta,'metrics100Pa',metrics);end
 clear mesh
end
report=array2table(rows,'VariableNames',{'level','max_height_displacement_m','max_slope_change_rad','max_curvature_change_per_m','relative_curvature_field_change'});writetable(report,fullfile(out,'pressure_convergence.csv'));
% Independent geometric check: exact parabola before/after a quadratic bend.
[XX,EE]=meshgrid(xi,eta);test=[XX(:)*.247,EE(:),.2*(XX(:)*.247).^2];du=zeros(size(test));du(:,3)=.1*test(:,1).^2;check=cad.posterior_metrics(test,du,xi,eta);x=.247*check.centers;exact=.6./(1+(.6*x).^2).^(3/2)-.4./(1+(.4*x).^2).^(3/2);error=max(abs(check.curvatureChange_per_m-exact),[],'all');assert(error<1e-9,'cad:Curvature','Curvature units/fit failed analytical parabola.');
saved=load(basisFile,'basis');B=kron(sample.H,speye(3))*saved.basis;d=load(fullfile(dynamicFolder,'U_10/dynamic.mat'),'dynamics','uvEquilibrium');times=[0,.0125,.025,.05];fields=cell(numel(times),1);curves=fields;
for k=1:numel(times),[~,index]=min(abs(d.dynamics.time-times(k)));q=d.uvEquilibrium.q+d.dynamics.Q(:,index);U=reshape(B*q,3,[])';fields{k}=struct('time',d.dynamics.time(index),'displacement',U);curves{k}=cad.posterior_metrics(sample.points,U,xi,eta);writematrix([sample.points,U],fullfile(out,"shape_frame_"+k+".csv"));end
scope=struct('xi',xi,'eta',eta,'pressureControl_Pa',100,'fitHalfWindow_over_c',.025,'curvatureDefinition',sample.metrics100Pa.definition,'analyticalParabolaMaximumError_per_m',error,'dynamicSource',dynamicFolder,'dynamicFrames','same original solid basis and UVLM trim as the fine time response; total displacement, not an imposed deformation','curvatureConvergence','reported separately; a fitted curvature is not a certified local stress or a pointwise edge derivative');cad.write_json(fullfile(out,'scope.json'),scope);save(fullfile(out,'fields.mat'),'sample','B','fields','curves','report','scope','-v7.3');
fig=figure('Visible','off','Color','w','Position',[100,100,1200,700]);if isprop(fig,'Theme'),fig.Theme='light';end;tiledlayout(2,1);nexttile;hold on;middle=ceil(numel(eta)/2);for k=1:numel(fields),U=reshape(fields{k}.displacement,numel(eta),numel(xi),3);plot(xi,1e3*U(middle,:,3),'DisplayName',sprintf('t=%.4g s',fields{k}.time));end;grid on;legend;xlabel('x/c material');ylabel('w_z [mm], total');title('Forma posterior del CAD, estación media; equilibrio + perturbación calculada');nexttile;hold on;for k=1:numel(curves),plot(curves{k}.centers,curves{k}.curvatureChange_per_m(middle,:),'DisplayName',sprintf('t=%.4g s',fields{k}.time));end;grid on;xlabel('x/c material');ylabel('cambio de curvatura [1/m]');title('Ajuste cúbico en ventana física; proyección cuerda/altura');exportgraphics(fig,fullfile(out,'posterior_shape_curvature.png'),'Resolution',150);close(fig);disp(report);
end
