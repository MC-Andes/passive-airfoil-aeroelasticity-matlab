function export_portal_fields(output)
%EXPORT_PORTAL_FIELDS Sample computed CAD motion for the MC-Andes website.
% Requires the completed local campaigns; performs no new physical simulation.
root=fileparts(fileparts(mfilename('fullpath')));
if ~isfolder(output),mkdir(output);end
p=jsondecode(fileread(fullfile(root,'results/cad_original/final_inputs.json')));
source=fullfile(root,'results/cad_original/posterior/fields.mat');
s=load(source,'sample','B');dynamic=fullfile(p.timeHalf128,'U_10/dynamic.mat');
if ~isfile(dynamic),dynamic=fullfile(root,dynamic);end
d=load(dynamic,'dynamics','uvEquilibrium');
indices=round(linspace(1,numel(d.dynamics.time),241));
ny=numel(s.sample.eta);middle=ceil(ny/2):ny:size(s.sample.points,1);Bz=s.B(3*middle,:);
writematrix([d.dynamics.time(indices)',1e6*(Bz*d.dynamics.Q(:,indices))'],fullfile(output,'posterior-motion.csv'));
writematrix([s.sample.xi(:),s.sample.points(middle,[1,3]),1e6*Bz*d.uvEquilibrium.q],fullfile(output,'posterior-reference.csv'));
scope=struct('sourceFields','results/cad_original/posterior/fields.mat','sourceFieldsSha256',cad.sha256(source), ...
    'sourceDynamics',erase(dynamic,[root,filesep]),'sourceDynamicsSha256',cad.sha256(dynamic), ...
    'columnsMotion','physical time [s], upper-skin displacement relative to UVLM trim [micrometres] at each xi', ...
    'columnsReference','xi=x/c, material x [m], material z [m], UVLM trim displacement [micrometres]', ...
    'samples',numel(indices),'eta',s.sample.eta(ceil(ny/2)),'scope','241 actual sampled states; no interpolated or imposed oscillation');
cad.write_json(fullfile(output,'motion-source.json'),scope);
end
