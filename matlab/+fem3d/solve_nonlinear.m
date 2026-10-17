function result = solve_nonlinear(system,load,bc,settings)
%SOLVE_NONLINEAR Continuation/Newton on the own finite-rotation solid.
% load may be a constant force vector, or @(u) -> [F,dFdu]. A follower-load
% callback MUST supply its consistent tangent. No undeclared dead/follower mix.
if nargin<4,settings=struct();end
o=struct('loadSteps',5,'maximumIterations',30,'relativeTolerance',1e-7,'forceScale',1e-8,'maximumLineSearch',15);
names=fieldnames(settings);for k=1:numel(names),assert(isfield(o,names{k}),'fem3d:Options','Unknown nonlinear setting.');o.(names{k})=settings.(names{k});end
nd=size(system.K,1);bc=fem3d.prepare_constraints(nd,bc);u=zeros(nd,1);history=[];
if ~isa(load,'function_handle'),if size(load,2)==3,load=reshape(load',[],1);end,load=load(:);end
for step=1:o.loadSteps
    factor=step/o.loadSteps; u(bc.dofs)=factor*bc.values;converged=false;
    for iteration=1:o.maximumIterations
        state=fem3d.nonlinear_assemble(system,u); [F,KL]=external(u);r=state.force-factor*F;
        metric=norm(r(bc.free),inf)/max(norm(factor*F(bc.free),inf),o.forceScale);
        history(end+1,:)=[step,iteration,metric,state.strainEnergy,state.minimumDeformationJacobian]; %#ok<AGROW>
        if metric<o.relativeTolerance,converged=true;break;end
        K=state.tangent-factor*KL; correction=zeros(nd,1);correction(bc.free)=-K(bc.free,bc.free)\r(bc.free);accepted=false;
        for search=0:o.maximumLineSearch-1
            candidate=u+2^(-search)*correction;
            try
                trial=fem3d.nonlinear_assemble(system,candidate);[nextF,~]=external(candidate);nextR=trial.force-factor*nextF;
                if norm(nextR(bc.free),inf)<(1-1e-4*2^(-search))*norm(r(bc.free),inf),u=candidate;accepted=true;break;end
            catch err
                if ~strcmp(err.identifier,'fem3d:Inversion'),rethrow(err);end
            end
        end
        assert(accepted,'fem3d:Nonlinear','Nonlinear line search failed at load factor %g.',factor);
    end
    assert(converged,'fem3d:Nonlinear','Nonlinear equilibrium failed at load factor %g.',factor);
end
result=struct('u',u,'displacement',reshape(u,3,[])','bc',bc,'history',history,'converged',converged,'state',state,'load',F);
    function [F,K]=external(displacement)
        if isa(load,'function_handle'),[F,K]=load(displacement);else,F=load;K=sparse(nd,nd);end
        assert(numel(F)==nd && isequal(size(K),[nd,nd]),'fem3d:Load','Nonlinear load/tangent dimensions invalid.');
    end
end
