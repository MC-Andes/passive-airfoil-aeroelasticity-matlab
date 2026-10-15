% INSPECT_REFERENCE Show converted CAD and extracted outer-envelope sections.
projectRoot=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(projectRoot,'matlab'));
folder=fullfile(projectRoot,'data','reference');
info=jsondecode(fileread(fullfile(folder,'cad_info.json')));
points=readmatrix(fullfile(folder,'cad_nodes.csv'));
triangles=readmatrix(fullfile(folder,'cad_triangles.csv'));
[surface,~]=pae.surface_from_sections(folder,4,10,"constant");
origin=surface.canonicalOriginCAD;
points=[points(:,1)-origin(1),points(:,3)-origin(2),points(:,2)-origin(3)];
triangles=triangles(:,[1,3,2]); % axis permutation reverses orientation
fig=figure('Color','w','Position',[80,80,1200,550]); tiledlayout(1,2);
if isprop(fig,'Theme'), fig.Theme='light'; end
nexttile; patch('Vertices',points,'Faces',triangles,'FaceColor',[.6,.7,.8], ...
    'EdgeColor','none','FaceAlpha',1); axis equal; view(35,25); camlight('headlight'); lighting gouraud;
xlim([min(points(:,1)),max(points(:,1))]); ylim([min(points(:,2)),max(points(:,2))]);
zlim([min(points(:,3)),max(points(:,3))]);
xlabel('Cuerda x [m]'); ylabel('Envergadura y [m]'); zlabel('Altura z [m]'); title('STEP: superficie triangulada');
nexttile; hold on;
for j=1:numel(info.sections)
    section=readmatrix(fullfile(folder,info.sections(j).file));
    plot(section(:,2),section(:,3),'DisplayName',sprintf('s=%.3f m',info.sections(j).span_m));
    plot(section(:,2),section(:,4),'HandleVisibility','off');
end
axis equal; grid on; xlabel('CAD X [m]'); ylabel('CAD Y [m]');
title('Envolventes exteriores: revisar antes de usarlas'); legend('Location','best');
axesHandles=findall(fig,'Type','axes');
set(axesHandles,'Color','w','XColor',[.15,.15,.15],'YColor',[.15,.15,.15], ...
    'ZColor',[.15,.15,.15],'GridColor',[.65,.65,.65],'FontSize',10);
outputFolder=fullfile(projectRoot,'results','reference'); if ~isfolder(outputFolder), mkdir(outputFolder); end
exportgraphics(fig,fullfile(outputFolder,'cad_inspection.png'),'Resolution',160);
disp(info);
