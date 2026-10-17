function x = factor_solve(L,order,b)
%FACTOR_SOLVE Apply an explicitly permuted sparse Cholesky inverse.
x=zeros(size(b)); x(order,:)=L'\(L\b(order,:));
end
