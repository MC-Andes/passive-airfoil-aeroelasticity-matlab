function snapshot = plot_circulatory_wake(tangentPath,dynamicPath,output)
%PLOT_CIRCULATORY_WAKE Reconstruct material emitted rows from retained states.
% Actual frozen-convection UVLM rings and their computed circulation; no
% invented velocity-deficit field. Prehistory far beyond the plot is omitted.
s=load(tangentPath,'fluid','model');d=load(dynamicPath,'dynamics','uvEquilibrium','flow');f=s.fluid;m=s.model;time=d.dynamics.time;[~,index]=min(abs(time-.1));dt=f.dt;g=f.state0.geometry;ambient=uvlm.air_velocity(d.flow,0,[0,0,0]);
if isfield(d.dynamics,'gammaTEPerturbation'),GTE=d.dynamics.gammaTEPerturbation;boundColor=zeros(g.nx*g.ny,1);else,GTE=f.TEextract*d.dynamics.gammaPerturbation;boundColor=d.dynamics.gammaPerturbation(:,index);end
maxAge=min(index-1,ceil(2*d.flow.referenceChord/(d.flow.speed*dt)));ny=f.spanPanels;X0=reshape(m.X0,[],3);indices=reshape(1:size(X0,1),g.nx+1,g.ny+1);te=indices(end,:);baseX=reshape(f.X0,[],3);wake=zeros(maxAge+1,ny+1,3);gamma=zeros(maxAge,ny);
for age=0:maxAge
 k=index-age;motion=reshape(m.PhiSurface*d.dynamics.Q(:,k),3,[])';wake(age+1,:,:)=reshape(baseX(te,:)+motion(te,:)+age*dt*ambient,1,ny+1,3);
 if age<maxAge
  emitted=max(1,k-1);gamma(age+1,:)=GTE(:,emitted)';
 end
end
[vertices,rings]=uvlm.wake_lattice(wake);X=reshape(baseX+reshape(m.PhiSurface*d.dynamics.Q(:,index),3,[])',size(f.X0));bound=uvlm.lattice(X);fig=figure('Visible','off','Color','w','Position',[100,100,1200,700]);if isprop(fig,'Theme'),fig.Theme='light';end
patch('Faces',bound.ring,'Vertices',bound.ringVertices,'FaceVertexCData',boundColor,'FaceColor',[.7,.7,.7],'EdgeColor',[.3,.3,.3]);hold on;
patch('Faces',rings,'Vertices',vertices,'FaceVertexCData',gamma(:),'FaceColor','flat','EdgeColor','none','FaceAlpha',.8);axis equal;grid on;view(35,25);xlabel('x [m]');ylabel('envergadura [m]');zlabel('altura [m], escala real');cb=colorbar;cb.Label.String='perturbación de circulación [m^2/s]';title(sprintf('CAD, U=%g m/s, t=%.4g s: anillos UVLM emitidos y convectados',d.flow.speed,time(index)));
if ~isfolder(output),mkdir(output);end
exportgraphics(fig,fullfile(output,'circulatory_wake.png'),'Resolution',160);close(fig);
snapshot=struct('time_s',time(index),'wakeRowsPlotted',maxAge,'historyRowsRetained',f.wakeRows,'convection','ambient frozen; actual TE material motion retained','verticalExaggeration',1,'field','circulation perturbation in m2/s, not viscous deficit','kelvin','every emitted row uses the previous bound trailing-edge circulation');cad.write_json(fullfile(output,'wake_scope.json'),snapshot);save(fullfile(output,'wake_snapshot.mat'),'wake','gamma','vertices','rings','snapshot');
end
