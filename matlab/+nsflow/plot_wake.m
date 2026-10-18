function plot_wake(result,output,titleText)
%PLOT_WAKE Actual computed mean velocity, deficit and vorticity on a declared plane.
if nargin<3,titleText='Navier–Stokes laminar propio: plano medio declarado';end
if isfolder(output),output=fullfile(output,'computed_wake.png');end
g=result.grid;d=g.dimension;U=result.meanVelocity;center=cell(1,d);body=[];if isfield(result,'solidMask'),body=result.solidMask;end
if isfield(result,'summary'),titleText=sprintf('%s; Re=%g, t*=[%g,%g]',titleText,result.summary.Re,result.summary.statisticsWindow);end
for i=1:d,center{i}=(U{i}+nsflow.shift(U{i},i,-1))/2;end
if d==2
    x=g.pressureCoordinates{1};z=g.pressureCoordinates{2};u=center{1};w=center{2};
    vort=(nsflow.shift(w,1,-1)-nsflow.shift(w,1,1))/(2*g.h(1))-(nsflow.shift(u,2,-1)-nsflow.shift(u,2,1))/(2*g.h(2));
else
    j=round((g.count(3)+1)/2);x=squeeze(g.pressureCoordinates{1}(:,:,j));z=squeeze(g.pressureCoordinates{2}(:,:,j));u=squeeze(center{1}(:,:,j));w=squeeze(center{2}(:,:,j));
    if ~isempty(body),body=squeeze(body(:,:,j));end
    vort=(circshift(w,-1,1)-circshift(w,1,1))/(2*g.h(1))-(circshift(u,-1,2)-circshift(u,1,2))/(2*g.h(2));
end
if ~isempty(body),u(body>.5)=NaN;vort(body>.5)=NaN;end
fig=figure('Visible','off','Color','w','Position',[100,100,1300,800]);if isprop(fig,'Theme'),fig.Theme='light';end
tiledlayout(fig,2,1,'Padding','compact');ax=nexttile;contourf(ax,x,z,1-u,24,'LineStyle','none');if ~isempty(body),hold(ax,'on');contour(ax,x,z,body,[.5,.5],'k','LineWidth',1);end;axis(ax,'equal','tight');colorbar(ax);xlabel(ax,'x/c');ylabel(ax,'altura/c');title(ax,'Déficit medio 1-u/U, calculado');
ax=nexttile;contourf(ax,x,z,vort,32,'LineStyle','none');if ~isempty(body),hold(ax,'on');contour(ax,x,z,body,[.5,.5],'k','LineWidth',1);end;axis(ax,'equal','tight');colorbar(ax);xlabel(ax,'x/c');ylabel(ax,'altura/c');title(ax,'Vorticidad transversal calculada, escala U/c');sgtitle(fig,titleText);
exportgraphics(fig,output,'Resolution',150);close(fig);
end
