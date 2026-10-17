function info = provenance(root,config,output)
%PROVENANCE Record exact inputs and code, including uncommitted work.
[status,commit]=system(['git -C ',cad.shell_quote(root),' rev-parse HEAD']); if status~=0, commit='unavailable'; end
[~,dirty]=system(['git -C ',cad.shell_quote(root),' status --short']);
info=struct('createdUTC',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')), ...
    'matlabVersion',version,'matlabRelease',version('-release'),'computer',computer, ...
    'commit',strtrim(commit),'gitStatus',dirty,'case',config,'transform', ...
    cad.transform(config.geometry.axisOrder(:)',config.geometry.origin_m(:)'));
paths=[dir(fullfile(root,'matlab','**','*.m'));dir(fullfile(root,'tools','*.py'));dir(fullfile(root,'tools','native','*.c'));dir(fullfile(root,'cases','**','*.json'))];
native=which('aero_vortex_kernel');if ~isempty(native),info.nativeKernel=struct('path',native,'sha256',cad.sha256(native),'scope','optional identical real Biot-Savart kernel; MATLAB fallback available');end
files=strings(numel(paths),1); hashes=files;
snapshot=fullfile(output,'executed_sources');if ~isfolder(snapshot),mkdir(snapshot);end
for k=1:numel(paths)
    p=fullfile(paths(k).folder,paths(k).name); files(k)=string(erase(p,[root,filesep])); hashes(k)=cad.sha256(p);
    destination=fullfile(snapshot,files(k));parent=fileparts(destination);if ~isfolder(parent),mkdir(parent);end;copyfile(p,destination);
end
info.sourceSnapshot='executed_sources: exact source content at campaign start; later source edits do not alter it';
writetable(table(files,hashes),fullfile(output,'source_hashes.csv'));
cad.write_json(fullfile(output,'provenance.json'),info);
end
