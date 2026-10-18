function [bc,info] = boundary_conditions(mesh,c,command)
%BOUNDARY_CONDITIONS Bonded shaft + ideal rigid driven front face, no root clamp.
assert(strcmp(mesh.info.source_sha256,c.geometry.sourceSha256),'morphing:Hash','Refusing labels from another STEP.');
labels=cad.face_labels(mesh);fixed=[];center=c.audit.shaft_socket.canonical_center_m;
for face=c.supports.cadFaces(:)'
 nodes=unique(mesh.faces(labels==face,:));assert(~isempty(nodes),'morphing:Support','Missing socket face %d.',face);
 X=mesh.nodes(nodes,:);r=X(:,[1,3])-center([1,3])';half=c.audit.shaft_socket.side_m/2;
 assert(all(abs(abs(r(:,1))-half)<1e-8 | abs(abs(r(:,2))-half)<1e-8),'morphing:Support','Socket signature mismatch.');
 fixed=[fixed;nodes(X(:,2)<=c.supports.insertion_m+1e-10 & X(:,2)>=-1e-10)]; %#ok<AGROW>
end
fixed=unique(fixed);drive=unique(mesh.faces(ismember(labels,c.slider.cadFaces),:));
assert(~isempty(fixed) && ~isempty(drive) && isempty(intersect(fixed,drive)),'morphing:BC','Missing or intersecting physical interfaces.');
assert(max(abs(mesh.nodes(drive,1)-(c.audit.slider_drive.native_plane_x_m-c.geometry.origin_m(1))))<1e-8,'morphing:Slider','Drive-face plane mismatch.');
dofs=[fem3d.dofs(fixed);fem3d.dofs(drive)];values=zeros(numel(dofs),1);values(3*numel(fixed)+(1:3:3*numel(drive)))=command;
bc=fem3d.prepare_constraints(3*size(mesh.nodes,1),struct('dofs',dofs,'values',values));fem3d.check_supports(mesh,bc);
info=struct('mountNodes',fixed,'driveNodes',drive,'insertionRequested_m',c.supports.insertion_m, ...
 'lastConstrainedSpan_m',max(mesh.nodes(fixed,2)),'command_m',command,'mountFaces',c.supports.cadFaces,'driveFaces',c.slider.cadFaces, ...
 'scope','three translations per solid node; rigid rotations restrained by face displacement field; all other root skin/tip nodes remain free unless on an actual selected interface');
end
