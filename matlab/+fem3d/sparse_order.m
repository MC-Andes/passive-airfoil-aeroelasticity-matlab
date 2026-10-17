function [order,info] = sparse_order(A)
%SPARSE_ORDER Choose sparse ordering by measured symbolic fill, not DOF count.
% MATLAB base AMD and nested dissection; no external structural solver.
p=symamd(A); first=sum(symbfact(A(p,p))); order=p; second=Inf;
if exist('dissect','file') || exist('dissect','builtin')
    candidate=dissect(A); second=sum(symbfact(A(candidate,candidate)));
    if second<first,order=candidate;end
end
info=struct('amdFactorNonzeros',first,'nestedDissectionFactorNonzeros',second, ...
    'predictedFactorNonzeros',min(first,second));
end
