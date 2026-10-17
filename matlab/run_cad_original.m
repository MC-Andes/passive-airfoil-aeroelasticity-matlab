function study = run_cad_original(casePath,runTet10)
%RUN_CAD_ORIGINAL Reproducible prescribed-load CAD structural campaign.
% Requires the three configured CSV meshes. Does not overwrite historical runs.
% MATLAB base only. Meshes are CAD preprocessing; all FEM is implemented here.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'matlab'));
if nargin<1 || isempty(casePath), casePath=fullfile(root,'cases','cad_original','case.json'); end
config=cad.load_case(root,casePath);
baseOrder=4; if isfield(config.discretization,'baseElementOrder'), baseOrder=config.discretization.baseElementOrder; end
assert(ismember(baseOrder,[4,10]),'cad:Case','Require element order 4 or 10.');
if nargin<2, runTet10=baseOrder==4; end
assert(~runTet10 || baseOrder==4,'cad:Case','The optional Tet10 pilot requires a Tet4 base campaign.');
label=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS'));
output=fullfile(root,'results','cad_original',['campaign_',label]);
assert(~isfolder(output),'cad:Output','Campaign output already exists.'); mkdir(output);
cad.provenance(root,config,output);
levels=string(config.discretization.meshes); count=numel(levels);
assert(count>=config.acceptance.minimumMeshLevels,'cad:Convergence','Need at least three meshes.');
for k=1:count
    assert(isfile(fullfile(root,levels(k),'mesh_info.json')),'cad:MissingMesh', ...
        'Prepare %s first; see docs/cad-original.md.',levels(k));
end
rows=cell(count+double(runTet10),1); stats=struct([]); observations=cell(size(rows)); modesCell=cell(size(rows));
for k=1:numel(rows)
    timer=tic;
    if k<=count, folder=levels(k); name="tet"+baseOrder+"_level_"+k; order=baseOrder;
    else, folder=string(config.discretization.tet10PilotMesh); name="tet10_pilot"; order=10; end
    out=fullfile(output,name); mkdir(out);
    fprintf('CAD %s: reading %s\n',name,folder);
    mesh=fem3d.read_mesh(fullfile(root,folder),config.geometry.axisOrder,config.geometry.origin_m);
    assert(strcmp(mesh.info.source_sha256,config.geometry.sourceSha256),'cad:Hash','Mesh source mismatch.');
    inputs=["nodes.csv","tetra.csv","faces.csv","mesh_info.json"]; hashes=strings(4,1);
    for j=1:4, hashes(j)=cad.sha256(fullfile(root,folder,inputs(j))); end
    writetable(table(inputs',hashes,'VariableNames',{'file','sha256'}),fullfile(out,'input_hashes.csv'));
    cad.atlas(mesh,config,fullfile(out,'atlas'));
    budget=cad.resource_estimate(mesh);
    assert(budget.assemblyChunkTripletBytes<=config.discretization.maxMatrixTripletBytes, ...
        'cad:Budget','Configured triplet budget is below the assembler chunk requirement.');
    cad.write_json(fullfile(out,'resources_before.json'),budget);
    if order==10
        assert(budget.tet10Dofs<=config.discretization.maximumPilotDofs,'cad:Budget', ...
            'Tet10 pilot exceeds configured DOF budget; measure resources before increasing it.');
        mesh=fem3d.elevate(mesh);
    end
    sensor=cad.sensors(mesh,config); patch=cad.pressure_patch(mesh,config);
    material=fem3d.material_reference(string(config.material.reference));
    fprintf('CAD %s: assemble %d DOFs, %d elements\n',name,3*size(mesh.nodes,1),size(mesh.tetra,1));
    ta=tic; system=fem3d.assemble(mesh,material); assemblySeconds=toc(ta);
    bc=struct('dofs',fem3d.dofs(find(mesh.nodes(:,1)<=config.supports.x_over_c*config.geometry.chord_m)), ...
        'description','ASSUMED anterior clamp; monolithic scenario A');
    ts=tic; result=fem3d.solve_static(system,patch.unitForce,bc); staticSeconds=toc(ts);
    observed=sensor.O*result.u; observations{k}=observed;
    bvec=reshape(patch.unitForce',[],1); scaleForce=max(norm(sum(patch.unitForce,1)),eps);
    forceRel=norm(result.forceBalance)/scaleForce;
    momentRel=norm(result.momentBalance)/(scaleForce*config.geometry.chord_m);
    unitForceRel=patch.forceDefect_N_per_Pa/scaleForce;
    unitMomentRel=patch.momentDefect_Nm_per_Pa/(scaleForce*config.geometry.chord_m);
    testU=mesh.nodes*[.1,.2,0;0,.1,0;.1,-.2,.3]'+[.1,.2,.3]; P=-patch.weights.*patch.normals;
    work1=sum(testU.*patch.unitForce,'all'); work2=sum((patch.H*testU).*P,'all');
    workRel=abs(work1-work2)/max([abs(work1),abs(work2),1e-12]);
    assert(max([unitForceRel,unitMomentRel,workRel])<=config.acceptance.transferRelative,'cad:Transfer','Transfer criterion failed.');
    fprintf('CAD %s: static residual %.3g; modes\n',name,result.relativeResidual);
    tm=tic; modes=fem3d.modal(system,bc,config.discretization.modalCount); modalSeconds=toc(tm);
    Oxyz=kron(sensor.H,speye(3)); modes.observed=Oxyz*modes.Phi; modesCell{k}=rmfield(modes,'Phi');
    modalTail=sensor.O*modes.Phi; tailRatio=zeros(numel(modes.frequencyHz),1);
    for j=1:numel(tailRatio)
        uj=reshape(modes.Phi(:,j),3,[])'; tailRatio(j)=max(abs(modalTail(:,j)))/max(vecnorm(uj,2,2));
    end
    writetable(table((1:numel(tailRatio))',modes.frequencyHz,modes.eigenResidual,tailRatio, ...
        'VariableNames',{'mode','frequency_Hz','eigen_residual','tail_z_to_max_amplitude'}),fullfile(out,'modes.csv'));
    % Same load distribution is scaled exactly in this linear model.
    pressures=config.loads.pressure_Pa(:); Y=observed*pressures';
    tableSensor=array2table([sensor.eta,sensor.points,observed],'VariableNames',{'eta','x_m','y_m','z_m','w_m_per_Pa'});
    writetable(tableSensor,fullfile(out,'sensors.csv'));
    writetable(array2table([pressures,result.maximumDisplacement*pressures, ...
        max(abs(observed))*pressures,result.strainEnergy*pressures.^2], ...
        'VariableNames',{'pressure_Pa','maximum_u_m','maximum_abs_sensor_w_m','strain_energy_J'}),fullfile(out,'pressure_response.csv'));
    reductionRows=[]; modalLoad=modes.Phi'*bvec;
    for n=unique([1,2,4,numel(modalLoad)])
        if n>numel(modalLoad), continue; end
        reduced=modalTail(:,1:n)*(modalLoad(1:n)./modes.omega(1:n).^2);
        reductionRows(end+1,:)=[n,norm(reduced-observed,inf)/max(norm(observed,inf),1e-15)]; %#ok<AGROW>
    end
    writetable(array2table(reductionRows,'VariableNames',{'retained_modes','static_sensor_relative_error'}),fullfile(out,'reduction.csv'));
    % Decoupled FRF is conditional on assumed modal damping and mode truncation.
    model=struct('K',diag(modes.omega.^2),'M',eye(numel(modalLoad)), ...
        'C',diag(2*config.discretization.dampingRatio*modes.omega),'B',modalLoad,'O',modalTail);
    frequency=linspace(0,1.25*modes.frequencyHz(end),401); frf=ael.harmonic(model,frequency,1);
    frfTable=array2table([frequency',abs(frf.Y)',rad2deg(angle(frf.Y))']);
    writetable(frfTable,fullfile(out,'conditional_modal_frf.csv'));
    cad.plot_response(sensor,pressures,Y,frf,config,out);
    rplot=result; rplot.displacement=result.displacement*max(pressures); rplot.stress.vonMises=result.stress.vonMises*max(pressures);
    fig=fem3d.plot_result(system,rplot,"CAD original | presión impuesta | apoyo/material supuestos | sin validación",20);
    exportgraphics(fig,fullfile(out,'deformation.png'),'Resolution',150); close(fig);
    % Compare the historical 16-neighbor footprint on identical point forces.
    if order==4
        regionNodes=unique(mesh.faces(ismember(cad.face_labels(mesh),config.regions.upperTailCadFaces),:));
        [Hnear,infoNear]=fem3d.transfer(mesh.nodes,patch.points,regionNodes,16);
        nearResult=fem3d.solve_static(system,Hnear'*P,bc);
        nearObserved=sensor.O*nearResult.u;
        transferComparison=struct('neighbors',16,'maximumRadius_m',infoNear.maximumPatchRadius, ...
            'meanRadius_m',mean(infoNear.patchRadii),'relativeSensorDifference', ...
            norm(nearObserved-observed,inf)/max(norm(observed,inf),1e-15), ...
            'interpretation','same integration point forces, same skin, different transfer footprint');
        cad.write_json(fullfile(out,'footprint_comparison.json'),transferComparison);
        clear nearResult Hnear;
    end
    mem=whos('system');
    stats(k).name=name; stats(k).order=order; stats(k).nodes=size(mesh.nodes,1); stats(k).tetra=size(mesh.tetra,1);
    stats(k).mass_kg=system.totalMass; stats(k).volumeError=mesh.info.relative_volume_error;
    stats(k).energy_J_per_Pa2=result.strainEnergy; stats(k).maxSensor_m_per_Pa=max(abs(observed));
    stats(k).residual=result.relativeResidual; stats(k).energyResidual=result.energyResidual; stats(k).maximumEigenEnergyResidual=max(modes.eigenEnergyResidual); stats(k).maximumEigenResidual=max(modes.eigenResidual); stats(k).forceBalanceRelative=forceRel; stats(k).momentBalanceRelative=momentRel;
    stats(k).area_m2=patch.area_m2; stats(k).projectedArea_m2=patch.projectedArea_m2;
    stats(k).assemblySeconds=assemblySeconds; stats(k).staticSeconds=staticSeconds; stats(k).modalSeconds=modalSeconds;
    stats(k).systemBytes=mem.bytes; stats(k).factorNonzeros=result.factorNonzeros; stats(k).totalSeconds=toc(timer);
    patch=rmfield(patch,{'H','unitForce'}); %#ok<NASGU>
    save(fullfile(out,'solution.mat'),'config','result','sensor','patch','modes','-v7.3');
    cad.write_json(fullfile(out,'summary.json'),stats(k));
    fprintf('CAD %s complete in %.1f s; sensor %.6g m/Pa\n',name,stats(k).totalSeconds,stats(k).maxSensor_m_per_Pa);
    clear system result mesh modes;
end
summary=struct2table(stats); writetable(summary,fullfile(output,'summary.csv'));
changes=struct([]);
for k=2:count
    changes(k-1).from=stats(k-1).name; changes(k-1).to=stats(k).name;
    changes(k-1).sensorRelative=norm(observations{k}-observations{k-1},inf)/max(norm(observations{k},inf),1e-15);
    changes(k-1).energyRelative=abs(stats(k-1).energy_J_per_Pa2/stats(k).energy_J_per_Pa2-1);
    changes(k-1).modes=cad.modal_compare(modesCell{k-1},modesCell{k},config.sensors.eta(:));
end
last=changes(end); relevant=true(size(last.modes.MAC));
if isfield(config.acceptance,'frequencyBand_Hz')
    freq=modesCell{count}.frequencyHz(last.modes.pairing);band=config.acceptance.frequencyBand_Hz;relevant=freq>=band(1) & freq<=band(2);
end
passed=last.sensorRelative<config.acceptance.tailEnergyRelative && ...
    last.energyRelative<config.acceptance.tailEnergyRelative && ...
    any(relevant) && max(last.modes.frequencyRelativeChange(relevant))<config.acceptance.frequencyRelative && ...
    all(last.modes.MAC(relevant)>.9) && ~any(last.modes.closeModes(relevant));
checks=all(summary.volumeError<config.acceptance.volumeRelative) && all(summary.energyResidual<config.acceptance.residualRelative) && ...
    all(summary.maximumEigenEnergyResidual<config.acceptance.residualRelative) && ...
    all(summary.forceBalanceRelative<config.acceptance.residualRelative) && ...
    all(summary.momentBalanceRelative<config.acceptance.residualRelative);
study=struct('output',output,'summary',stats,'meshChanges',changes,'elementOrder',baseOrder,'discretizationCriteriaMet',passed, ...
    'algebraicCriteriaMet',checks,'structuralModelClosed',passed && checks,'experimentalValidation',false,'relevantModeMask',relevant,'numericalClosure',passed && checks, ... ...
    'scope','Conditional monolithic CAD prescribed-load study; numerical convergence and physical validation assessed separately.', ...
    'blockedStages','Physical mount/material/slider require evidence; two-way FSI/UVLM/dynamic and viscous wake gates remain open.');
if runTet10
    study.elementComparison=struct('sensorRelative',norm(observations{end}-observations{1},inf)/max(norm(observations{end},inf),1e-15), ...
        'energyRelative',abs(stats(1).energy_J_per_Pa2/stats(end).energy_J_per_Pa2-1), ...
        'geometry','same faceted CAD mesh, straight-edge elevation');
end
cad.write_json(fullfile(output,'acceptance.json'),study);
fprintf('CAD campaign: %s\nDiscretization convergence criteria met: %d. Conditional numerical model closed: %d.\n',output,passed,passed && checks);
end
