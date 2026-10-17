function plane = plot_cross_wake(result,output,xChord)
%PLOT_CROSS_WAKE Computed finite-span crossflow and streamwise vorticity.
if nargin<3,xChord=1.25;end
g=result.grid;assert(g.dimension==3,'nsflow:Plane','A three-dimensional field is required.');
[~,ix]=min(abs(g.pressureCoordinates{1}(:,1,1)-xChord));xActual=g.pressureCoordinates{1}(ix,1,1);
U=cell(1,3);for k=1:3,U{k}=(result.meanVelocity{k}+nsflow.shift(result.meanVelocity{k},k,-1))/2;end
height=squeeze(g.pressureCoordinates{2}(ix,:,:));span=squeeze(g.pressureCoordinates{3}(ix,:,:));
deficit=1-squeeze(U{1}(ix,:,:));v=squeeze(U{2}(ix,:,:));w=squeeze(U{3}(ix,:,:));
curl=(nsflow.shift(U{3},2,-1)-nsflow.shift(U{3},2,1))/(2*g.h(2))-(nsflow.shift(U{2},3,-1)-nsflow.shift(U{2},3,1))/(2*g.h(3));curl=squeeze(curl(ix,:,:));
fig=figure('Visible','off','Color','w','Position',[100,100,1200,650]);if isprop(fig,'Theme'),fig.Theme='light';end
layout=tiledlayout(fig,1,2,'Padding','loose','TileSpacing','loose');layout.OuterPosition=[0,0,1,.86];nexttile;contourf(span,height,deficit,24,'LineStyle','none');hold on;
step=4;take={1:step:size(v,1),1:step:size(v,2)};quiver(span(take{:}),height(take{:}),w(take{:}),v(take{:}),.8,'k');axis equal tight;colorbar;xlabel('envergadura/c');ylabel('altura/c');title({'Déficit 1-u/U y flujo transversal','vectores con autoescala'},'FontSize',12);
nexttile;contourf(span,height,curl,24,'LineStyle','none');axis equal tight;colorbar;xlabel('envergadura/c');ylabel('altura/c');title({'Vorticidad en dirección de corriente','escala U/c'},'FontSize',12);
titleText=sprintf('CAD 3D laminar; x/c=%.4g; campos medios calculados',xActual);if isfield(result,'summary'),titleText=sprintf('%s; Re=%g; t*=[%g,%g]',titleText,result.summary.Re,result.summary.statisticsWindow);end
annotation(fig,'textbox',[.01,.91,.98,.07],'String',titleText,'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
if ~isfolder(output),mkdir(output);end;exportgraphics(fig,fullfile(output,'cross_wake.png'),'Resolution',150);close(fig);
plane=struct('requestedXOverChord',xChord,'actualXOverChord',xActual,'span',span,'height',height,'deficit',deficit,'spanVelocity',w,'heightVelocity',v,'streamwiseVorticity',curl,'scope','computed finite-time laminar mean; capped exterior; no PIV or turbulent wake claim');save(fullfile(output,'cross_wake.mat'),'plane');
end
