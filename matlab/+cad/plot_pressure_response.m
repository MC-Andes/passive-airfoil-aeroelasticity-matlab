function report = plot_pressure_response(tangentPath,dynamicPath,output)
%PLOT_PRESSURE_RESPONSE Equivalent thin-surface mean load and initial inertia.
% The Bernoulli field is only the computed initial gamma-dot contribution,
% not the complete fluctuating skin Cp and not a prescribed sinusoid.
s=load(tangentPath,'fluid','model');d=load(dynamicPath,'dynamics','uvEquilibrium','flow');f=s.fluid;m=s.model;g=f.state0.geometry;X=reshape(f.X0,[],3);
meanPressure=f.base.pressureJump;gammaDot=f.E*d.dynamics.A(:,1);inertialPressure=d.flow.rho*gammaDot;
nodal=g.Hcent'*(inertialPressure.*g.panels.area.*g.panels.normal);observed=m.PhiSurface'*reshape(nodal',[],1);expected=f.Ld*gammaDot;
defect=norm(observed-expected)/max(norm(expected),realmin);assert(defect<1e-11,'cad:Pressure','Initial Bernoulli pressure projection is inconsistent.');
span=max(X(:,2))-min(X(:,2));vertices=[X(:,1)/d.flow.referenceChord,(X(:,2)-min(X(:,2)))/span,zeros(size(X,1),1)];faces=g.panels.corners(:,[1,2,4,3]);
fig=figure('Visible','off','Color','w','Position',[100,100,1200,620]);if isprop(fig,'Theme'),fig.Theme='light';end
layout=tiledlayout(fig,1,2,'Padding','loose','TileSpacing','loose');fields={meanPressure,inertialPressure};labels={{'Carga media equivalente normal','no es Cp de ambas pieles'},{'Contribución Bernoulli inicial','rho * dGamma/dt, por aceleración calculada'}};
for k=1:2,nexttile;patch('Faces',faces,'Vertices',vertices,'FaceVertexCData',fields{k},'FaceColor','flat','EdgeColor','none');view(2);axis tight;xlabel('x/c');ylabel('eta envergadura');cb=colorbar;cb.Label.String='Pa';title(labels{k},'FontSize',12);end
title(layout,sprintf('CAD UVLM, U=%g m/s: presiones equivalentes calculadas',d.flow.speed),'FontSize',13);
if ~isfolder(output),mkdir(output);end;exportgraphics(fig,fullfile(output,'equivalent_pressure.png'),'Resolution',150);close(fig);
points=g.Hcent*X;writetable(array2table([points,meanPressure,inertialPressure,g.panels.area],'VariableNames',{'x_m','span_m','height_m','mean_equivalent_pressure_Pa','initial_Bernoulli_pressure_Pa','area_m2'}),fullfile(output,'pressure_fields.csv'));
report=struct('maximumMeanAbsolute_Pa',max(abs(meanPressure)),'maximumInitialBernoulliAbsolute_Pa',max(abs(inertialPressure)), ...
 'initialBernoulliProjectionRelativeDefect',defect,'UVLMEquilibriumCL',d.uvEquilibrium.aero.CL,'UVLMEquilibriumMaximumTE_m',max(abs(m.observation*d.uvEquilibrium.q)), ...
 'scope','mean equivalent normal resultant per panel; initial acceleration-induced Bernoulli contribution only; actual Joukowski force arms retained by solver; no separately resolved skin Cp');cad.write_json(fullfile(output,'scope.json'),report);
end
