function value = divergence(g,U)
%DIVERGENCE Finite-volume face flux divergence at cell centers.
value=zeros(g.count);
for i=1:g.dimension,value=value+(nsflow.shift(U{i},i,-1)-U{i})/g.h(i);end
end
