function r = solve_contact(system,F,bc,normal,gap,C,d)
%SOLVE_CONTACT Frictionless affine unilateral contact by an exact active set.
% normal*u+gap >= 0, contact reaction normal'*lambda, lambda >= 0.
% Small-motion gap model only; interfaces must be supplied, never inferred.
nd=size(system.K,1);if nargin<6,C=sparse(0,nd);d=zeros(0,1);end
gap=gap(:);nc=size(normal,1);active=false(nc,1);
assert(size(normal,2)==nd && numel(gap)==nc,'fem3d:Contact','Invalid contact gap operator.');
for iteration=1:100
    eq=[C;normal(active,:)];target=[d(:);-gap(active)];
    solution=fem3d.solve_constrained(system,F,bc,eq,target);gaps=normal*solution.u+gap;
    lambda=zeros(nc,1);lambda(active)=-solution.lambda(size(C,1)+(1:sum(active)));
    release=find(active & lambda<-1e-9,1);close=find(~active & gaps<-1e-10,1);
    if ~isempty(release),active(release)=false;
    elseif ~isempty(close),active(close)=true;
    else
        r=solution;r.gap=gaps;r.normalReaction=lambda;r.active=active;r.contactIterations=iteration;
        r.complementarity=max(abs(lambda.*gaps));return;
    end
end
error('fem3d:Contact','Contact active set did not converge.');
end
