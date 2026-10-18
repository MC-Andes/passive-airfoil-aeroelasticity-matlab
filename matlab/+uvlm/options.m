function o = options(user)
%OPTIONS Explicit UVLM numerical choices; no viscous or turbulence model.
o=struct('coreRadius',1e-5,'wakeLength',20,'convection','frozen', ...
    'convectionOrder',2,'courantLimit',1,'maximumWakeRows',5000,'steadyStart',true);
if nargin>0
    names=fieldnames(user);
    for k=1:numel(names)
        assert(isfield(o,names{k}),'uvlm:Options','Unknown option %s.',names{k}); o.(names{k})=user.(names{k});
    end
end
assert(o.coreRadius>0 && o.wakeLength>0 && ismember(string(o.convection),["frozen","free"]) && ...
    ismember(o.convectionOrder,[1,2]) && o.maximumWakeRows>=2,'uvlm:Options','Invalid vortex/wake options.');
end
