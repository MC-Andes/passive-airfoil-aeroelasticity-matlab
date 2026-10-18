function result = uvlm_static(nodes,velocity,flow,settings)
%UVLM_STATIC Stationary ring/wake limit, with explicit matching dynamic options.
if nargin<4,settings=struct('wakeLength',1000*flow.referenceChord);end
[~,result]=uvlm.initialize(nodes,velocity,flow,settings);
end
