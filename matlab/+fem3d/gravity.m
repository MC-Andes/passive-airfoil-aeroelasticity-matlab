function F = gravity(system,acceleration)
%GRAVITY Consistent body load from the same mass operator used in dynamics.
assert(numel(acceleration)==3 && all(isfinite(acceleration)),'fem3d:Gravity','Finite Cartesian acceleration required.');
F=reshape(system.M*repmat(acceleration(:),size(system.mesh.nodes,1),1),3,[])';
end
