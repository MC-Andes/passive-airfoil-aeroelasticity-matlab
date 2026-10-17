function value = sha256(path)
%SHA256 Byte-exact fingerprint using OS utilities; no JVM required.
% Hashing only: no numerical computation is delegated to shell utilities.
assert(isfile(path),'cad:File','Cannot read %s.',path);
if ispc
    assert(~contains(path,'"'),'cad:Path','Invalid Windows filename.');
    [status,out]=system(['certutil -hashfile "',char(path),'" SHA256']);
    match=regexp(out,'(?m)^\s*([0-9a-fA-F]{64})\s*$','tokens','once');
else
    quote=char(39); safe=[quote,strrep(char(path),quote,[quote,'"',quote,'"',quote]),quote];
    [status,out]=system(['shasum -a 256 ',safe]);
    if status~=0, [status,out]=system(['sha256sum ',safe]); end
    match=regexp(out,'^([0-9a-fA-F]{64})\s','tokens','once');
end
assert(status==0 && ~isempty(match),'cad:Hash','SHA256 utility unavailable.'); value=lower(match{1});
end
