function result = linearize(model,solution,overrides)
%LINEARIZE Local QUASI-STEADY dynamics about an equilibrium.
% Eigenvalues are diagnostics of the instantaneous-load model. They do not
% establish a flutter boundary: wake memory and added mass are absent.
if nargin<3, overrides=struct(); end
settings=pae.options(overrides);
assert(solution.converged,'pae:Linearization','A converged equilibrium is required.');
q=solution.q; r=model.ndof; flow=solution.state.flow;
Qq=pae.jacobian(@(x) load_at(x,zeros(r,1)),q,settings.differenceStep);
Qv=pae.jacobian(@(v) load_at(q,v),zeros(r,1),settings.differenceStep);
[~,K]=pae.internal_force(model,q);
effectiveK=K-Qq; effectiveC=model.structure.C-Qv;
A=[zeros(r),eye(r);-model.structure.M\effectiveK,-model.structure.M\effectiveC];
lambda=eig(A);
result=struct('A',A,'eigenvalues',lambda,'Qaero_q',Qq,'Qaero_qdot',Qv, ...
    'effectiveK',effectiveK,'effectiveC',effectiveC,'maximumRealPart',max(real(lambda)), ...
    'scope','Local quasi-steady diagnostic. Not a validated flutter prediction.');
    function Q=load_at(x,v)
        state=pae.evaluate(model,x,v,flow); Q=state.Qaero;
    end
end
