function V = air_velocity(flow,t,P)
%AIR_VELOCITY Inertial ambient velocity, including optional prescribed gust.
beta=0; if isfield(flow,'beta'), beta=flow.beta; end
V=repmat(flow.speed*[cos(flow.alpha)*cos(beta),sin(beta),sin(flow.alpha)*cos(beta)],size(P,1),1);
if isfield(flow,'gust'), V=V+flow.gust(t,P); end
assert(isequal(size(V),size(P)) && all(isfinite(V(:))),'uvlm:Flow','Invalid ambient velocity.');
end
