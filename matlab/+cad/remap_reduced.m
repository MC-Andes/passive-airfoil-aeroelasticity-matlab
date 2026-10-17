function bundle = remap_reduced(source,basisFile,nx,ny,output)
%REMAP_REDUCED Change only the fluid discretization; structural basis fixed.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));s=load(source,'bundle');bundle=s.bundle;c=bundle.config;
path=fullfile(root,string(c.discretization.meshes{end}));mesh=fem3d.read_mesh(path,c.geometry.axisOrder,c.geometry.origin_m);if size(mesh.tetra,2)==4,mesh=fem3d.elevate(mesh);end
p=load(basisFile,'basis');surface=cad.surface(mesh,c,nx,ny);bundle.surface=surface;bundle.model.X0=surface.nodes;bundle.model.PhiSurface=surface.Oxyz*p.basis;
if ~isfolder(output),mkdir(output);end
save(fullfile(output,'reduced_cad.mat'),'bundle','-v7.3');
end
