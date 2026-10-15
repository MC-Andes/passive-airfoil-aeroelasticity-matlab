function indices = dofs(nodes,components)
%DOFS One-based structural DOFs; components 1=x, 2=y, 3=z.
arguments
    nodes double {mustBeInteger,mustBePositive}
    components double {mustBeInteger,mustBePositive} = 1:3
end
assert(all(components(:)<=3),'fem3d:DOF','Components must be 1,2,3.');
indices=reshape((3*nodes(:)-3+components(:)')',[],1);
end
