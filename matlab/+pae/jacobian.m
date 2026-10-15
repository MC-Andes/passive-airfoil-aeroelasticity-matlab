function J = jacobian(fun,x,step)
%JACOBIAN Central finite differences; step is in generalized coordinate units.
x=x(:); baseline=fun(x); J=zeros(numel(baseline),numel(x));
if isscalar(step), step=repmat(step,numel(x),1); else, step=step(:); end
assert(numel(step)==numel(x) && all(step>0 & isfinite(step)), ...
    'pae:Difference','Invalid finite-difference steps.');
for j=1:numel(x)
    delta=zeros(size(x)); delta(j)=step(j)*max(1,abs(x(j)));
    J(:,j)=(fun(x+delta)-fun(x-delta))/(2*delta(j));
end
end
