function [x,w] = gauss01(n)
%GAUSS01 Gauss-Legendre nodes/weights on [0,1], without toolboxes.
assert(isscalar(n) && n==fix(n) && n>=1 && n<=12,'fem3d:Quadrature','Order must be 1..12.');
b=(1:n-1)./sqrt(4*(1:n-1).^2-1);
[V,D]=eig(diag(b,1)+diag(b,-1)); [x,i]=sort(diag(D));
w=V(1,i)'.^2; x=(x+1)/2;
end
