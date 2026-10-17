function [tri,owners] = plot_faces(faces)
%PLOT_FACES Subdivide Tri6 for display only; preserve physical face ownership.
nf=size(faces,1);
if size(faces,2)==3, tri=faces; owners=(1:nf)'; return; end
assert(size(faces,2)==6,'fem3d:Face','Expected Tri3/Tri6.');
tri=[faces(:,[1,4,6]);faces(:,[4,2,5]);faces(:,[6,5,3]);faces(:,[4,5,6])];
owners=repmat((1:nf)',4,1);
end
