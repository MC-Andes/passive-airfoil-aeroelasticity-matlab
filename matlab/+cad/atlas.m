function report = atlas(mesh,config,output)
%ATLAS Faceted geometry signatures, local material thickness and element locations.
if ~isfolder(output), mkdir(output); end
X=mesh.nodes; F=mesh.faces; labels=cad.face_labels(mesh); tags=unique(labels);
a=X(F(:,1),:); b=X(F(:,2),:); c=X(F(:,3),:); ar=.5*vecnorm(cross(b-a,c-a,2),2,2); ctr=(a+b+c)/3;
rows=zeros(numel(tags),12); group=strings(numel(tags),1);
for k=1:numel(tags)
    take=labels==tags(k); ids=unique(F(take,:)); lo=min(X(ids,:),[],1); hi=max(X(ids,:),[],1);
    rows(k,:)=[tags(k),sum(take),sum(ar(take)),sum(ar(take).*ctr(take,:),1)/sum(ar(take)),lo,hi];
    group(k)="unclassified boundary; mechanical role unknown";
    if ismember(tags(k),config.regions.upperTailCadFaces), group(k)="outer upper tail candidate; patch exposure tested"; end
    if ismember(tags(k),config.regions.lowerTailCadFaces), group(k)="outer lower tail candidate; no pressure assigned"; end
    if ismember(tags(k),config.regions.spanEndCadFaces), group(k)="span end cap"; end
end
signatures=array2table(rows,'VariableNames',{'cad_face','triangles','area_m2','cx_m','cy_m','cz_m', ...
    'xmin_m','ymin_m','zmin_m','xmax_m','ymax_m','zmax_m'}); signatures.group=group;
writetable(signatures,fullfile(output,'face_signatures.csv'));
% Quality is tied to physical locations, not only a global minimum.
T=mesh.tetra(:,1:4); ne=size(T,1); ratio=zeros(ne,1); vol=ratio; cent=zeros(ne,3);
for e=1:ne
    xe=X(T(e,:),:); J=(xe(2:4,:)-xe(1,:))'; sv=svd(J);
    ratio(e)=sv(end)/sv(1); vol(e)=abs(det(J))/6; cent(e,:)=mean(xe,1);
end
[~,ix]=sort(ratio); ix=ix(1:min(200,ne));
worst=table(ix,ratio(ix),vol(ix),cent(ix,1),cent(ix,2),cent(ix,3), ...
    'VariableNames',{'element','singular_ratio','volume_m3','x_m','y_m','z_m'});
writetable(worst,fullfile(output,'worst_elements.csv'));
% Vertical skin thickness is a directional material interval, not normal thickness.
xi=linspace(.35,.995,41); eta=config.sensors.eta(:); tab=zeros(numel(xi)*numel(eta),8); k=0;
for j=1:numel(eta), for x=xi
    k=k+1; h=cad.ray(mesh,[x*config.geometry.chord_m,(eta(j)-.5)*config.geometry.span_m]); z=h.z; n=numel(z);
    tab(k,:)=[x,eta(j),n,NaN,NaN,NaN,NaN,NaN];
    if n>=2 && rem(n,2)==0
        tab(k,4:8)=[z(end)-z(1),z(end)-z(end-1),z(2)-z(1),sum(z(2:2:end)-z(1:2:end)),1];
    end
end, end
thickness=array2table(tab,'VariableNames',{'x_over_c','eta','crossings','envelope_m', ...
    'upper_material_interval_z_m','lower_material_interval_z_m','total_solid_intervals_z_m','even_crossings'});
writetable(thickness,fullfile(output,'directional_thickness.csv'));
report=struct('faceCount',numel(tags),'minimumJacobianRatio',min(ratio),'elementsBelow1e4',sum(ratio<1e-4), ...
    'volume_m3',sum(vol),'thicknessInterpretation','faceted vertical material intervals; not local normal skin thickness', ...
    'unresolved','mount, rib/flexure mechanical roles, slider interfaces, normal thickness and experimental evidence');
cad.write_json(fullfile(output,'atlas.json'),report);
fig=figure('Visible','off','Color','w','Position',[100,100,1200,500]);
if isprop(fig,'Theme'), fig.Theme='light'; end
subplot(1,2,1); hold on;
for e=config.regions.upperTailCadFaces(:)'
    use=labels==e; patch('Faces',F(use,1:3),'Vertices',X,'FaceColor',[.3,.65,.85],'EdgeColor','none');
end
scatter3(cent(ix,1),cent(ix,2),cent(ix,3),9,log10(ratio(ix)),'filled');
axis equal; grid on; view(35,25); xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]'); title('Piel superior y 200 elementos de peor calidad'); cb=colorbar; cb.Label.String='log10(razón singular)';
subplot(1,2,2); mid=abs(tab(:,2)-.5)<1e-8;
plot(tab(mid,1),1e3*tab(mid,4:6),'LineWidth',1.2); grid on;
xlabel('x/c'); ylabel('Intervalo vertical [mm]'); legend('Envolvente','Piel superior','Piel inferior','Location','best');
title('Sección central: intervalos de material');
exportgraphics(fig,fullfile(output,'atlas.png'),'Resolution',150); close(fig);
end
