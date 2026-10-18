function out = shell_quote(value)
%SHELL_QUOTE Quote paths passed to fixed read-only provenance commands.
value=char(value);
if ispc
    assert(~any(ismember(value,['"','%','!'])),'cad:Path','Unsupported Windows shell path.'); out=['"',value,'"'];
else
    q=char(39); out=[q,strrep(value,q,[q,'"',q,'"',q]),q];
end
end
