function mesh = read_mesh(folder,axisOrder,origin)
%READ_MESH Read cached STEP volume mesh; optionally permute axes and translate.
% origin is expressed in the permuted axes, not in original STEP axes.
arguments
    folder
    axisOrder (1,3) double = [1,2,3]
    origin (1,3) double = [0,0,0]
end
assert(isequal(sort(axisOrder),1:3),'fem3d:Axes','Axis order must be a permutation.');
mesh.info=jsondecode(fileread(fullfile(folder,'mesh_info.json')));
assert(strcmp(mesh.info.units,'m'),'fem3d:Units','Mesh must be in metres.');
if isfield(mesh.info,'relative_volume_error')
    assert(mesh.info.relative_volume_error<.01,'fem3d:MeshVolume', ...
        'Mesh volume differs from the CAD audit by more than 1%%; inspect geometry/preprocessing.');
end
X=readmatrix(fullfile(folder,'nodes.csv')); T=readmatrix(fullfile(folder,'tetra.csv'));
F=readmatrix(fullfile(folder,'faces.csv'));
mesh.nodes=X(:,axisOrder)-origin; assert(ismember(size(T,2),[5,11]) && size(F,2)==(size(T,2)-1)/2+2, ...
    'fem3d:Mesh','Require Tet4/Tri3 or Tet10/Tri6 CSVs with CAD labels.');
mesh.tetra=T(:,1:end-1);
mesh.cadVolume=T(:,end); mesh.cadFaces=F(:,1:end-1); mesh.cadFaceId=F(:,end);
% Compute orientations from volume connectivity; CAD face orientation may differ.
[mesh.faces,mesh.faceOwners]=fem3d.boundary(mesh.tetra,mesh.nodes);
mesh.materialId=ones(size(T,1),1);
mesh.source="STEP volume mesh (preprocessing only), "+string(mesh.info.source_name);
end
