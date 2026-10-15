function [Q,nodalLoads] = map_loads(J,Hforce,panelForces)
%MAP_LOADS Adjoint load map preserves virtual work exactly for this model.
nodalLoads=Hforce'*panelForces;
Q=J'*nodalLoads(:);
end
