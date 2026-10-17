function result = vlm(nodes,velocity,flow)
%VLM Conservative nodal adapter for the stationary horseshoe solver.
result=pae.vlm(nodes,velocity,flow);
result.nodalForces=result.panels.Hforce'*result.panelForces;
end
