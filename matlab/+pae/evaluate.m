function state = evaluate(model,q,qdot,flow)
%EVALUATE Structural force and aerodynamics on the current geometry.
q=q(:); qdot=qdot(:);
assert(numel(qdot)==model.ndof && all(isfinite(qdot)),'pae:State','Invalid velocity.');
[nodes,J]=pae.geometry(model,q);
velocity=reshape(J*qdot,[],3);
if ~isfield(flow,'referenceArea'), flow.referenceArea=model.surface.referenceArea; end
if ~isfield(flow,'referenceChord'), flow.referenceChord=model.surface.referenceChord; end
if isfield(model.structure,'aeroFcn')
    aero=model.structure.aeroFcn(nodes,velocity,flow);
else
    aero=pae.vlm(nodes,velocity,flow,model.aeroOptions);
end
assert(isfield(aero,'panelForces') && isequal(size(aero.panelForces),[model.surface.nx*model.surface.ny,3]) ...
    && all(isfinite(aero.panelForces(:))),'pae:Aero','Invalid aerodynamic loads.');
panel=pae.panels(nodes);
[Q,nodalLoads]=pae.map_loads(J,panel.Hforce,aero.panelForces);
[force,tangent]=pae.internal_force(model,q);
state=struct('nodes',nodes,'J',J,'aero',aero,'Qaero',Q,'nodalLoads',nodalLoads, ...
    'internalForce',force,'internalTangent',tangent,'flow',flow);
end
