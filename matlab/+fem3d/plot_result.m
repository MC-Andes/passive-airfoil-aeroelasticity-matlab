function fig = plot_result(system,result,titleText,scale)
%PLOT_RESULT Supports, actual displacement magnitude and element stresses.
arguments
    system struct
    result struct
    titleText (1,1) string
    scale (1,1) double {mustBePositive} = 1
end
X=system.mesh.nodes; [F,faceMap]=fem3d.plot_faces(system.mesh.faces); U=result.displacement;
fig=figure('Color','w','WindowStyle','normal','Position',[100,100,1500,700],'Visible','off');
if isprop(fig,'Theme'), fig.Theme='light'; end
annotation(fig,'textbox',[.04,.92,.92,.06],'String',titleText, ...
    'Interpreter','none','FontSize',12,'HorizontalAlignment','center','EdgeColor','none');
annotation(fig,'textbox',[.05,.80,.4,.05], ...
    'String',sprintf('Desplazamiento; dibujo amplificado x%.0f',scale), ...
    'FontSize',11,'HorizontalAlignment','center','EdgeColor','none');
annotation(fig,'textbox',[.55,.80,.4,.05],'String','Tensión por elemento en la superficie', ...
    'FontSize',11,'HorizontalAlignment','center','EdgeColor','none');
ax=axes(fig,'Position',[.06,.20,.36,.56]); hold(ax,'on');
patch(ax,'Faces',F,'Vertices',X+scale*U,'FaceVertexCData',1e3*vecnorm(U,2,2), ...
    'FaceColor','interp','EdgeColor','none');
fixed=unique(ceil(result.bc.dofs/3)); fixed=fixed(1:max(1,ceil(numel(fixed)/350)):end);
scatter3(ax,X(fixed,1),X(fixed,2),X(fixed,3),7,'k','filled');
axis(ax,'equal'); grid(ax,'on'); view(ax,35,25); xlabel(ax,'x cuerda [m]');
ylabel(ax,'y envergadura [m]'); zlabel(ax,'z altura [m]');
cb=colorbar(ax); cb.Label.String='|u| [mm]'; ax.FontSize=10;
cb.Position=[.45,.25,.015,.46];
ax=axes(fig,'Position',[.56,.20,.36,.56]); vals=result.stress.vonMises(system.mesh.faceOwners(faceMap))/1e6;
patch(ax,'Faces',F,'Vertices',X,'FaceVertexCData',vals,'FaceColor','flat','EdgeColor','none');
axis(ax,'equal'); grid(ax,'on'); view(ax,35,25); xlabel(ax,'x cuerda [m]');
ylabel(ax,'y envergadura [m]'); zlabel(ax,'z altura [m]');
cb=colorbar(ax); cb.Label.String='von Mises [MPa]'; ax.FontSize=10;
cb.Position=[.95,.25,.015,.46];
end
