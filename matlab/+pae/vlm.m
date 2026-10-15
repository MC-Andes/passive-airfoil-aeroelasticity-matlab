function result = vlm(nodes,nodeVelocity,flow,options)
%VLM Steady/quasi-steady 3D horseshoe vortex-lattice aerodynamic loads.
% No viscous drag, separation, aerodynamic memory, or added-mass force.
% Positive alpha means a positive upward freestream component in body axes.
% panelForces act at the quarter-chord bound-vortex midpoint. pressureJump
% is an equivalent thin-surface load, NOT resolved skin pressure.
arguments
    nodes double
    nodeVelocity double
    flow struct
    options struct = struct()
end
p = pae.panels(nodes); np = numel(p.area);
assert(isequal(size(nodeVelocity),[numel(nodes)/3,3]) && all(isfinite(nodeVelocity(:))), ...
    'pae:Velocity','Invalid aerodynamic surface velocities.');
assert(isfield(flow,'rho') && isfield(flow,'speed') && isfield(flow,'alpha'), ...
    'pae:Flow','Flow requires rho, speed and alpha (radians).');
assert(isscalar(flow.speed) && flow.speed>=0 && isfinite(flow.speed) && ...
    isscalar(flow.rho) && flow.rho>0 && isfinite(flow.rho) && ...
    isscalar(flow.alpha) && isfinite(flow.alpha),'pae:Flow','Invalid flow.');
beta = 0; if isfield(flow,'beta'), beta = flow.beta; end
assert(isscalar(beta) && isfinite(beta),'pae:Flow','Invalid sideslip.');
direction = [cos(flow.alpha)*cos(beta),sin(beta),sin(flow.alpha)*cos(beta)];
ref = [0,0,0]; if isfield(flow,'referencePoint'), ref = flow.referencePoint(:)'; end
assert(numel(ref)==3 && all(isfinite(ref)),'pae:Flow','Invalid reference point.');
scale = max(vecnorm(p.B-p.A,2,2));
extent = max(max(reshape(nodes,[],3))-min(reshape(nodes,[],3)));
core = get_option(options,'coreRadius',max(extent,scale)*1e-6);
wakeLength = get_option(options,'wakeLength',100*extent);
assert(core>0 && wakeLength>0 && isfinite(core+wakeLength),'pae:VLM','Invalid vortex options.');
forces = zeros(np,3); gamma = zeros(np,1); residual = 0; condition = 1;
if flow.speed>0
    influence = zeros(np,np);
    velocityAtForce = zeros(np,3,np);
    for j=1:np
        farA = p.A(j,:)+wakeLength*direction;
        farB = p.B(j,:)+wakeLength*direction;
        vc = pae.segment_velocity(p.control,farA,p.A(j,:),core) ...
           + pae.segment_velocity(p.control,p.A(j,:),p.B(j,:),core) ...
           + pae.segment_velocity(p.control,p.B(j,:),farB,core);
        influence(:,j) = sum(vc.*p.normal,2);
        vf = pae.segment_velocity(p.forcePoint,farA,p.A(j,:),core) ...
           + pae.segment_velocity(p.forcePoint,p.A(j,:),p.B(j,:),core) ...
           + pae.segment_velocity(p.forcePoint,p.B(j,:),farB,core);
        velocityAtForce(:,:,j) = vf;
    end
    condition = rcond(influence);
    assert(condition>1e-12,'pae:VLM','Ill-conditioned aerodynamic lattice; inspect geometry/mesh.');
    relative = flow.speed*direction-p.Hcontrol*nodeVelocity;
    rhs = -sum(relative.*p.normal,2);
    gamma = influence\rhs;
    residual = norm(influence*gamma-rhs,inf);
    induced = zeros(np,3);
    for j=1:np, induced = induced+velocityAtForce(:,:,j)*gamma(j); end
    relativeForce = flow.speed*direction-p.Hforce*nodeVelocity+induced;
    forces = flow.rho*gamma.*cross(relativeForce,p.B-p.A,2);
end
total = sum(forces,1);
moment = sum(cross(p.forcePoint-ref,forces,2),1);
liftDirection = cross(direction,[0,1,0]); liftDirection = liftDirection/norm(liftDirection);
lift = dot(total,liftDirection); drag = dot(total,direction);
area = get_option(flow,'referenceArea',sum(p.area));
chord = get_option(flow,'referenceChord',extent);
assert(isscalar(area) && area>0 && isfinite(area) && isscalar(chord) && chord>0 && isfinite(chord), ...
    'pae:Flow','Invalid coefficient reference dimensions.');
qS = .5*flow.rho*flow.speed^2*area;
CL=NaN; CDi=NaN; Cm=NaN;
if qS>0, CL=lift/qS; CDi=drag/qS; Cm=moment(2)/(qS*chord); end
result = struct('panelForces',forces,'gamma',gamma,'totalForce',total, ...
    'totalMoment',moment,'lift',lift,'inducedDrag',drag,'CL',CL,'CDi',CDi,'Cm',Cm, ...
    'pressureJump',sum(forces.*p.normal,2)./p.area, ...
    'boundaryResidual',residual,'influenceRcond',condition,'panels',p, ...
    'fidelity','quasi-steady inviscid horseshoe VLM; not UVLM');
end

function value=get_option(s,name,fallback)
value=fallback; if isfield(s,name), value=s.(name); end
end
