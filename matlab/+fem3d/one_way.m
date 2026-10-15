function result = one_way(system,surface,flow,bc,candidates)
%ONE_WAY Decoupled steady aerodynamics -> conservative forces -> linear FEM.
% Aerodynamics evaluated ONCE on initial shape, no aerodynamic feedback.
arguments
    system struct
    surface struct
    flow struct
    bc struct
    candidates double = unique(system.mesh.faces(:))
end
flow.referenceArea=surface.referenceArea; flow.referenceChord=surface.referenceChord;
nv=numel(surface.nodes)/3;
aero=pae.vlm(surface.nodes,zeros(nv,3),flow);
[H,transferInfo]=fem3d.transfer(system.mesh.nodes,aero.panels.forcePoint,candidates);
F=H'*aero.panelForces;
structure=fem3d.solve_static(system,F,bc);
forceDefect=sum(F,1)-aero.totalForce;
momentDefect=sum(cross(system.mesh.nodes,F,2),1)- ...
    sum(cross(aero.panels.forcePoint,aero.panelForces,2),1);
pointDisplacement=H*structure.displacement;
virtualWorkDefect=sum(pointDisplacement.*aero.panelForces,'all')-structure.u'*reshape(F',[],1);
result=struct('aero',aero,'structure',structure,'H',H,'transferInfo',transferInfo, ...
    'transferForceDefect',forceDefect,'transferMomentDefect',momentDefect, ...
    'virtualWorkDefect',virtualWorkDefect,'flow',flow, ...
    'coupling','one-way on initial geometry; linear structural solid FEM');
end
