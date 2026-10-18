function r = solve_friction(system,F,bc,normal,gap,tangent,mu,C,d,previous)
%SOLVE_FRICTION One normal and one tangential guide coordinate, Coulomb law.
% Quasi-static rate-independent increment; previous.tangentialPosition stores slip.
% General 3D contact surfaces and finite sliding need geometric updates.
nd=size(system.K,1);if nargin<9,C=sparse(0,nd);d=zeros(0,1);end
if nargin<10,previous=struct('tangentialPosition',0);end
xold=previous.tangentialPosition;
assert(size(normal,1)==1 && size(tangent,1)==1 && mu>=0,'fem3d:Friction','One guide contact and nonnegative friction required.');
if size(F,2)==3,F=reshape(F',[],1);else,F=F(:);end
stick=fem3d.solve_contact(system,F,bc,normal,gap,[C;tangent],[d(:);xold]);
normalForce=stick.normalReaction;required=-stick.lambda(size(C,1)+1);
if abs(required)<=mu*normalForce+1e-9
    r=stick;r.state='stick';r.frictionForce=required;r.dissipatedWork=0;r.tangentialPosition=tangent*r.u;return;
end
free=fem3d.solve_contact(system,F,bc,normal,gap,C,d);direction=sign(tangent*free.u-xold);if direction==0,direction=-sign(required);end
force=-mu*free.normalReaction*direction;
for iteration=1:100
    solution=fem3d.solve_contact(system,F+tangent'*force,bc,normal,gap,C,d);
    target=-mu*solution.normalReaction*direction;
    if abs(target-force)<1e-8*max(1,abs(target)),break;end,force=(force+target)/2;
end
assert(iteration<100,'fem3d:Friction','Coulomb guide iteration failed.');
r=solution;r.state='slip';r.frictionForce=target;r.tangentialPosition=tangent*solution.u;r.dissipatedWork=-target*(r.tangentialPosition-xold);
assert(r.dissipatedWork>=-1e-12,'fem3d:Friction','Friction supplied energy instead of dissipating it.');
end
