function result = run_tail_response()
%RUN_TAIL_RESPONSE Trailing-edge motion under mean VLM + prescribed oscillation.
% Flat plate surrogate, NOT the CAD/slider. Harmonic pressure is synthetic.
% The steady VLM does not generate pressure fluctuations or a physical wake.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'matlab'));
out=fullfile(root,'results','tail_response'); if ~isfolder(out), mkdir(out); end
c=.247; b=.36; xf=.3*c;
mat=struct('E',2.11e9,'nu',.3,'rho',1180); % nominal published PLA model
x=unique([linspace(0,c,33),xf]); y=linspace(-b/2,b/2,33);
mesh=plate.mesh(x,y,@(x,y) .004-.002*(x>xf)); %#ok<INUSD>
s=plate.assemble(mesh,mat);
fixed=fem3d.dofs(find(mesh.nodes(:,1)<=xf+1e-12));
span=linspace(-b/2,b/2,21)'; points=[c*ones(size(span)),span];
O=plate.probe(mesh,points,1); Otheta=plate.probe(mesh,points,2);
plan=struct('spanStations',[-b/2;b/2],'leadingEdge',[0,-b/2,0;0,b/2,0],'chord',[c;c]);
surface=pae.make_surface(plan,8,16);
flow=struct('rho',1.225,'speed',10,'alpha',3*pi/180, ...
    'referenceArea',c*b,'referenceChord',c);
aero=pae.vlm(surface.nodes,zeros(numel(surface.nodes)/3,3),flow);
B=plate.pressure_map(mesh,surface.nodes(:,1,1),surface.nodes(1,:,2));
pmean=aero.pressureJump;
static=ael.static(s.K,fixed,B*pmean,O);
% Compare modal truncation and time step, without changing pressure field.
model=ael.reduce(s.K,s.M,fixed,24,.02,B,O);
modalStatic=ael.harmonic(model,0,pmean);
modalStaticError=norm(modalStatic.Y-static.Y)/max(norm(static.Y),eps);
excitationHz=.8*model.frequencyHz(1); pAmplitude=.1*pmean;
frequencies=linspace(.1*model.frequencyHz(1),1.5*model.frequencyHz(4),350);
frf=ael.harmonic(model,frequencies,pAmplitude);
duration=25/model.frequencyHz(1); time=linspace(0,duration,2501);
pressure=pmean+pAmplitude*sin(2*pi*excitationHz*time);
inputFcn=ael.pressure_history(time,pressure);
% Begin at modal static equilibrium; subsequent motion comes from imposed input.
q0=model.K\(model.B*pmean);
dynamic=ael.transient(model,time,inputFcn,q0);
coarse=ael.transient(model,time(1:2:end),inputFcn,q0);
timeStepError=norm(coarse.Y-dynamic.Y(:,1:2:end),'fro')/max(norm(dynamic.Y(:,1:2:end),'fro'),eps);
fewer=model; keep=1:12;
for field=["M","C","K"], fewer.(field)=model.(field)(keep,keep); end
fewer.B=model.B(keep,:); fewer.O=model.O(:,keep);
frf12=ael.harmonic(fewer,excitationHz,pAmplitude);
frf24=ael.harmonic(model,excitationHz,pAmplitude);
modeCountError=norm(frf12.Y-frf24.Y)/max(norm(frf24.Y),eps);
window=dynamic.time>=duration/2;
metrics=ael.response_metrics(dynamic.time(window),dynamic.Y(:,window));
rmsDisplacement=metrics.fluctuationRMS;
rotation=Otheta*static.U;
% Keep load interpretation explicit: uniform panel pressure acts at panel centroids.
F=B*pmean; distributedMoment=sum(mesh.nodes(:,1).*F(1:3:end));
summary=table(model.frequencyHz(1),excitationHz,max(abs(static.Y)),max(rmsDisplacement), ...
    max(abs(rotation)),modalStaticError,timeStepError,modeCountError, ...
    sum(F(1:3:end)),aero.totalForce(3),distributedMoment,-aero.totalMoment(2), ...
    'VariableNames',{'first_mode_Hz','imposed_excitation_Hz','max_TE_static_m', ...
    'max_TE_fluctuation_RMS_m','max_TE_rotation_rad','modal_static_relative_error', ...
    'time_step_relative_change','mode_count_relative_change','mapped_vertical_force_N', ...
    'VLM_vertical_force_N','distributed_pressure_chord_moment_Nm','VLM_chord_moment_Nm'});
writetable(summary,fullfile(out,'summary.csv'));
writetable(table(span,static.Y,rotation,rmsDisplacement,'VariableNames', ...
    {'span_m','TE_static_m','TE_rotation_rad','TE_fluctuation_RMS_m'}),fullfile(out,'trailing_edge.csv'));
center=ceil(numel(span)/2);
writetable(table(dynamic.time',dynamic.Y(center,:)','VariableNames', ...
    {'time_s','TE_center_m'}),fullfile(out,'time_history.csv'));
writetable(table(frequencies',abs(frf.Y(center,:))',angle(frf.Y(center,:))','VariableNames', ...
    {'frequency_Hz','TE_center_peak_amplitude_m','TE_center_phase_rad'}),fullfile(out,'forced_response.csv'));
status=struct('geometry','flat plate surrogate of reference dimensions; no slider', ...
    'pressure','steady VLM mean plus imposed 10 percent sinusoidal pressure', ...
    'damping','2 percent modal damping assumed, not measured', ...
    'coupling','one-way prescribed loads on undeformed geometry', ...
    'wake','not computed: no unsteady or viscous wake solver implemented', ...
    'cad_convergence','original tetra4 CAD calculation remains unconverged', ...
    'pressure_moment','uniform panel pressure and VLM quarter-chord forces have different moments');
fid=fopen(fullfile(out,'scope.json'),'w'); cleaner=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(status,'PrettyPrint',true));
result=struct('summary',summary,'status',status,'static',static,'dynamic',dynamic, ...
    'model',model,'system',s,'aero',aero,'pressure',pressure,'frf',frf,'metrics',metrics);
save(fullfile(out,'result.mat'),'result','-v7.3');
fig=newfigure(); plot(span*1000,static.Y*1000,'LineWidth',1.3); grid on;
xlabel('Envergadura (mm)'); ylabel('Desplazamiento TE (mm)'); title('Borde de salida: carga media');
exportgraphics(fig,fullfile(out,'01_trailing_edge.png'),'Resolution',160); close(fig);
fig=newfigure(); plot(dynamic.time,dynamic.Y(center,:)*1000); grid on;
xlabel('Tiempo (s)'); ylabel('Desplazamiento TE (mm)'); title('Presión oscilante impuesta');
exportgraphics(fig,fullfile(out,'02_vibration.png'),'Resolution',160); close(fig);
fig=newfigure(); plot(frequencies,abs(frf.Y(center,:))*1000); grid on;
xlabel('Frecuencia (Hz)'); ylabel('Amplitud TE (mm)'); title('Respuesta a presión armónica');
exportgraphics(fig,fullfile(out,'03_frequency_response.png'),'Resolution',160); close(fig);
disp(summary);
end

function fig=newfigure()
fig=figure('Color','w','Visible','off','WindowStyle','normal');
if isprop(fig,'Theme'), fig.Theme='light'; end
end
