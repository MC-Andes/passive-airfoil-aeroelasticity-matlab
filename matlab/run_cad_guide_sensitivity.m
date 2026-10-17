function report = run_cad_guide_sensitivity(campaign)
%RUN_CAD_GUIDE_SENSITIVITY Conditional preload/friction and position families.
root=fileparts(fileparts(mfilename('fullpath')));s=load(fullfile(campaign,'mechanisms/guide_compliance.mat'),'response','mechanism','c');r=s.response;C=r.observation;G=r.compliance;Y=r.influence;u0=r.base(:,1);sensor=load(fullfile(campaign,'tet10_level_3/solution.mat'),'sensor');O=sensor.sensor.O;output=fullfile(campaign,'mechanisms');rows=[];
for offset=[-.0001,0,.0001]
 d=zeros(6,1);d(1)=offset;lambda=G\(C*u0-d);u=u0-Y*lambda;rows(end+1,:)=[offset,max(abs(O*u)),norm(C*u-d,inf)]; %#ok<AGROW>
end
writematrix(rows,fullfile(output,'position_offsets.csv'));
springs=[];for stiffness=[0,1e3,1e4,1e5,1e6]
 Ks=diag([ones(1,3),ones(1,3)*s.c.geometry.chord_m^2])*stiffness;lambda=(eye(6)+Ks*G)\(Ks*(C*u0));u=u0-Y*lambda;springs(end+1,:)=[stiffness,max(abs(O*u)),norm(C*u)]; %#ok<AGROW>
end
writematrix(springs,fullfile(output,'guide_stiffness.csv'));
% Condensation is exact for the base linear FEM. Contact coordinates are an
% explicit idealized guide with one normal and one chord friction coordinate.
small=struct('K',G\eye(6));small.K=(small.K+small.K')/2;bc=struct('dofs',[2,4,5,6]);normal=sparse([0,0,1,0,0,0]);tangent=sparse([1,0,0,0,0,0]);cycles=[];
for mu=[0,.2,.4]
 previous=struct('tangentialPosition',0);dissipation=0;
 for scale=[linspace(0,1,21),linspace(1,-1,41),linspace(-1,0,21)]
  F=small.K*(C*u0)*scale;F(3)=F(3)-1;F(1)=F(1)+scale;contact=fem3d.solve_friction(small,F,bc,normal,0,tangent,mu,sparse(0,6),zeros(0,1),previous);dissipation=dissipation+contact.dissipatedWork;
  cycles(end+1,:)=[mu,scale,contact.tangentialPosition,contact.normalReaction,contact.frictionForce,dissipation];previous=contact; %#ok<AGROW>
 end
end
writematrix(cycles,fullfile(output,'conditional_friction_cycles.csv'));
report=struct('positionOffsets_m',[-.0001,0,.0001],'preload_N',1,'imposedTangentialCycleAmplitude_N',1,'frictionCoefficients',[0,.2,.4],'dissipationNonnegative',all(cycles(:,end)>=-1e-12),'hypothesis','small-motion exact FEM compliance, idealized mean-coordinate guide; no measured surface friction, preload or clearance','notMeasuredHysteresis',true);cad.write_json(fullfile(output,'guide_sensitivity.json'),report);
fig=figure('Visible','off','Color','w','Position',[100,100,1100,700]);if isprop(fig,'Theme'),fig.Theme='light';end;tiledlayout(2,1);nexttile;plot(1e3*rows(:,1),1e3*rows(:,2),'o-');xlabel('desplazamiento de guía impuesto [mm]');ylabel('máximo |w_TE| [mm]');grid on;nexttile;hold on;for mu=[0,.2,.4],take=cycles(:,1)==mu;plot(1e3*cycles(take,3),cycles(take,5),'DisplayName',sprintf('mu=%g',mu));end;xlabel('coordenada tangencial [mm]');ylabel('fricción [N]');legend;grid on;exportgraphics(fig,fullfile(output,'guide_sensitivity.png'),'Resolution',150);close(fig);
end
