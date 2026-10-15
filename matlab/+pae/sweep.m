function [summary,solutions] = sweep(model,speeds,angles,flowBase,overrides)
%SWEEP Continuation over speed and incidence. Angles are radians.
if nargin<5, overrides=struct(); end
assert(isvector(speeds) && all(isfinite(speeds)) && all(speeds>=0) && ...
    isvector(angles) && all(isfinite(angles)),'pae:Sweep','Invalid sweep axes.');
count=numel(speeds)*numel(angles); solutions=cell(count,1);
values=zeros(count,8+model.ndof); row=0;
for alpha=angles(:)'
    q=zeros(model.ndof,1);
    for speed=speeds(:)'
        row=row+1; flow=flowBase; flow.speed=speed; flow.alpha=alpha;
        solution=pae.solve_static(model,flow,q,overrides);
        solutions{row}=solution;
        if solution.converged, q=solution.q; end
        displacement=reshape(solution.state.nodes-model.surface.nodes,[],3);
        a=solution.state.aero;
        values(row,:)=[speed,rad2deg(alpha),a.CL,a.CDi,a.Cm, ...
            max(vecnorm(displacement,2,2)),solution.residualNorm,solution.converged,solution.q'];
    end
end
names=[{'U_m_s','alpha_deg','CL','CDi','Cm','maxDisplacement_m','residual','converged'}, ...
    cellstr(matlab.lang.makeUniqueStrings(matlab.lang.makeValidName("q_"+string(model.structure.labels(:)'))))];
summary=array2table(values,'VariableNames',names);
summary.structuralSource=repmat(string(model.structure.source),count,1);
summary.coordinateUnits=repmat(strjoin(string(model.structure.coordinateUnits),';'),count,1);
summary.aerodynamicModel=repmat(string(solutions{1}.state.aero.fidelity),count,1);
end
