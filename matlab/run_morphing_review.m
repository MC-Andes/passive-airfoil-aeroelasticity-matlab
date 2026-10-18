function review = run_morphing_review(manifestPath)
%RUN_MORPHING_REVIEW Explicit new-geometry evidence; no inherited acceptance.
root=fileparts(fileparts(mfilename('fullpath')));p=jsondecode(fileread(manifestPath));out=fullfile(root,'results/morphing_wing/review');if ~isfolder(out),mkdir(out);end
c=morphing.load_case(root);cad.provenance(root,c,out);a=jsondecode(fileread(fullfile(p.tet4,'summary.json')));b=jsondecode(fileread(fullfile(p.tet10,'summary.json')));n=jsondecode(fileread(fullfile(p.nonlinear,'summary.json')));v=jsondecode(fileread(p.verification));
assert(strcmp(a.geometrySha256,c.geometry.sourceSha256) && strcmp(b.geometrySha256,c.geometry.sourceSha256) && strcmp(n.geometrySha256,c.geometry.sourceSha256),'morphing:Hash','Cross-geometry evidence rejected.');
review=struct('geometrySha256',c.geometry.sourceSha256,'inputs',p,'verification',v,'complementaryAnalysesRetained',{c.analysis.complementary}, ...
 'computed',struct('tet4',a,'tet10',b,'smallNonlinearTet4',n),'comparison',struct('TE_relative',abs(a.maxTE_m/b.maxTE_m-1),'servoForceRatio',a.driveForce_N(1)/b.driveForce_N(1),'energyRatio',a.energy_J/b.energy_J,'f1Ratio',a.frequencyHz(1)/b.frequencyHz(1),'smallNonlinear_vs_Tet4Tangent_TE_relative',abs(n.maxTE_m/(a.maxTE_m*n.command_m/a.command_m)-1)), ...
 'mechanicalMeshConvergenceEstablished',false,'tunnelCFDExecuted',false,'flutterBoundaryEstablished',false,'strengthQualified',false, ...
 'scope','updated study definition and verified implementation; primary tunnel-Re CFD and mechanical mesh convergence remain unestablished; old CAD evidence retained separately');
fields={'tet4','tet10'};interfaces=struct();
for k=1:2
 s=load(fullfile(p.(fields{k}),'solution.mat'),'interface');i=s.interface;interfaces.(fields{k})=struct('mountNodeCount',numel(i.mountNodes),'driveNodeCount',numel(i.driveNodes),'insertionRequested_m',i.insertionRequested_m,'lastConstrainedSpan_m',i.lastConstrainedSpan_m);
end
review.interfaces=interfaces;cad.write_json(fullfile(out,'review.json'),review);copyfile(manifestPath,fullfile(out,'inputs.json'));
files={manifestPath,fullfile(p.tet4,'summary.json'),fullfile(p.tet10,'summary.json'),fullfile(p.nonlinear,'summary.json'),p.verification,c.geometry.audit,c.studyData,'cases/morphing_wing/case.json'};hashes=strings(numel(files),1);for k=1:numel(files),hashes(k)=cad.sha256(files{k});end;writetable(table(string(files(:)),hashes,'VariableNames',{'file','sha256'}),fullfile(out,'evidence_hashes.csv'));
fig=figure('Visible','off','Color','w','Position',[100,100,1100,650]);if isprop(fig,'Theme'),fig.Theme='light';end
tiledlayout(2,1);nexttile;plot(c.sensors.eta,1e6*a.sensor_m,'o-',c.sensors.eta,1e6*b.sensor_m,'x-');xlabel('η');ylabel('w TE [µm]');grid on;legend('Tet4','Tet10','Location','best');nexttile;bar([a.frequencyHz,b.frequencyHz]);xlabel('Modo por orden frecuencial, slider bloqueado');ylabel('f [Hz]');legend('Tet4','Tet10','Location','northwest');grid on;exportgraphics(fig,fullfile(out,'element_comparison.png'),'Resolution',150);close(fig);
fid=fopen(fullfile(out,'report.md'),'w');cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Geometría actualizada: revisión de los cambios\n\n');fprintf(fid,'El estudio principal usa PETG, eje adherido y slider impuesto. Se mantienen vibraciones, flutter y resistencia como complementos.\n\n');fprintf(fid,'[Definición y reproducción](%s) · [Resultados técnicos](%s) · [Evidencia y hashes](%s)\n\n',fullfile(root,'docs/morphing-wing.md'),fullfile(root,'docs/morphing-wing-results.md'),fullfile(out,'evidence_hashes.csv'));
fprintf(fid,'![Secciones reales del CAD actualizado; rojo alojamiento, azul frente accionado](%s)\n\n',fullfile(out,'updated_sections.png'));
fprintf(fid,'La geometría tiene 630 caras y no es prismática. Cuerda nominal 250 mm; extensión CAD aproximadamente 247 mm, sin reescalado.\n\n');fprintf(fid,'![Comparación de elementos para el mismo apoyo y comando de -1 µm](%s)\n\n',fullfile(out,'element_comparison.png'));
fprintf(fid,'Tet4/Tet10: cambio TE %.3f%%; razón de reacción %.3f, energía %.3f y primera frecuencia %.3f. No se acepta convergencia mecánica con un solo cambio pequeño de TE.\n\n',100*review.comparison.TE_relative,review.comparison.servoForceRatio,review.comparison.energyRatio,review.comparison.f1Ratio);
fprintf(fid,'Tet10: %.0f DOFs, masa %.6f kg, TE máximo %.6f µm a -1 µm impuesto, f1 %.5f Hz. Valores condicionales y no convergidos.\n\n',b.dofs,b.mass_kg,1e6*b.maxTE_m,b.frequencyHz(1));
fprintf(fid,'Prueba no lineal Tet4 a %.3f mm: convergió, mínimo det(F) %.6f, máximo de deformación Green %.6g y máximo von Mises Cauchy %.6g MPa. Es un control de implementación, no la predicción de resistencia de la pieza.\n\n',1e3*n.command_m,n.minimumDeformationJacobian,n.maxGreenStrain,n.maximumCauchyVonMises_Pa/1e6);
fprintf(fid,'Pasan 75 controles. La tabla de 231 configuraciones es planificación; no contiene polares CFD calculadas. No se ha establecido flutter ni resistencia experimental. Los recorridos -12,5/-25 mm siguen siendo referencias propuestas, no límites confirmados.\n');
disp(review.comparison);
end
