function solution = solve_static(model,flow,q0,overrides)
%SOLVE_STATIC Two-way static equilibrium with damped Newton iterations.
if nargin<3 || isempty(q0), q0=zeros(model.ndof,1); end
if nargin<4, overrides=struct(); end
settings=pae.options(overrides); q=q0(:);
history=zeros(settings.maxIterations,3); converged=false; message='Iteration limit reached.';
for iteration=1:settings.maxIterations
    state=pae.evaluate(model,q,zeros(model.ndof,1),flow);
    load=state.Qaero+model.structure.externalForce;
    residual=state.internalForce-load;
    tolerance=settings.absoluteTolerance+settings.relativeTolerance*max(norm(load),norm(state.internalForce));
    history(iteration,:)=[iteration,norm(residual),norm(q)];
    if norm(residual)<=tolerance, converged=true; message='Equilibrium converged.'; break; end
    aeroTangent=pae.jacobian(@aerodynamic_load,q,settings.differenceStep);
    K=state.internalTangent-aeroTangent;
    referenceScale=max(norm(state.internalTangent,2),norm(aeroTangent,2));
    singularValues=svd(K);
    if rcond(K)<1e-12 || singularValues(end)<=1e-10*referenceScale
        message='Singular equilibrium tangent; no unique Newton correction.'; break;
    end
    correction=-K\residual;
    accepted=false;
    for search=0:settings.maxLineSearch-1
        candidate=q+2^(-search)*correction;
        next=pae.evaluate(model,candidate,zeros(model.ndof,1),flow);
        nextResidual=next.internalForce-next.Qaero-model.structure.externalForce;
        if norm(nextResidual)<(1-1e-4*2^(-search))*norm(residual)
            q=candidate; accepted=true; break;
        end
    end
    if ~accepted, message='Newton line search failed.'; break; end
end
state=pae.evaluate(model,q,zeros(model.ndof,1),flow);
residual=state.internalForce-state.Qaero-model.structure.externalForce;
tolerance=settings.absoluteTolerance+settings.relativeTolerance*max(norm(state.internalForce),norm(state.Qaero+model.structure.externalForce));
if norm(residual)<=tolerance, converged=true; message='Equilibrium converged.'; end
solution=struct('q',q,'converged',converged,'message',message,'iterations',iteration, ...
    'residualNorm',norm(residual),'history',history(1:iteration,:), ...
    'state',state,'structuralSource',model.structure.source);
if ~converged && settings.errorOnFailure
    error('pae:Equilibrium','%s Residual = %.6g.',message,norm(residual));
end
    function force=aerodynamic_load(coordinates)
        evaluated=pae.evaluate(model,coordinates,zeros(model.ndof,1),flow);
        force=evaluated.Qaero;
    end
end
