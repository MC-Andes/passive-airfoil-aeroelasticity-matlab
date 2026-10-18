function result = static(model,flow,settings)
%STATIC Two-way aeroelastic equilibrium in a CAD-derived structural basis.
% model.K/M/C, PhiSurface (3Na x Nm), X0, and aerodynamic nodal force callback.
% Independent load/displacement residuals, relaxation and explicit failure.
if nargin<3,settings=struct();end
o=struct('maximumIterations',150,'relaxation',.5,'loadTolerance',1e-5, ...
    'displacementTolerance',1e-5,'initialQ',[],'forceScale',1e-6,'displacementScale',1e-9,'errorOnFailure',true);
names=fieldnames(settings); for k=1:numel(names),assert(isfield(o,names{k}),'fsi:Options','Unknown setting.');o.(names{k})=settings.(names{k});end
n=size(model.K,1); q=zeros(n,1); if ~isempty(o.initialQ),q=o.initialQ(:);end; lastLoad=zeros(n,1); history=zeros(o.maximumIterations,4); converged=false;
for k=1:o.maximumIterations
    X=reshape(reshape(model.X0,[],3)'+reshape(model.PhiSurface*q,3,[]),3,[])';
    X=reshape(X,size(model.X0));
    aero=model.aero(X,zeros(numel(X)/3,3),flow);
    force=model.PhiSurface'*reshape(aero.nodalForces',[],1);
    equilibrium=model.K\force; delta=equilibrium-q;
    ru=norm(delta)/max([norm(equilibrium),norm(q),o.displacementScale]);
    rf=norm(force-lastLoad)/max([norm(force),norm(lastLoad),o.forceScale]);
    residual=norm(model.K*q-force)/max(norm(force),o.forceScale);
    history(k,:)=[ru,rf,residual,norm(q)];
    if ru<o.displacementTolerance && rf<o.loadTolerance && residual<o.loadTolerance
        converged=true; break;
    end
    q=q+o.relaxation*delta; lastLoad=force;
end
result=struct('q',q,'aero',aero,'load',force,'history',history(1:k,:), ...
    'converged',converged,'iterations',k,'relativeEquilibriumResidual',residual,'settings',o);
if o.errorOnFailure,assert(converged,'fsi:Equilibrium','Static FSI failed; residual %.3g.',residual);end
end
