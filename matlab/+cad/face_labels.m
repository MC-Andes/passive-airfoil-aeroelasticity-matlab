function ids = face_labels(mesh)
%FACE_LABELS Join outward FEM faces to CAD labels by connectivity, never proximity.
assert(isfield(mesh,'cadFaces') && isfield(mesh,'cadFaceId'),'cad:Labels','CAD face labels required.');
[ok,loc]=ismember(sort(mesh.faces(:,1:3),2),sort(mesh.cadFaces(:,1:3),2),'rows');
assert(all(ok) && numel(unique(loc))==size(mesh.faces,1),'cad:Labels','Incomplete/duplicate CAD boundary labels.');
ids=mesh.cadFaceId(loc);
end
