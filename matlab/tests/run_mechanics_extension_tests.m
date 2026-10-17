function report = run_mechanics_extension_tests()
%RUN_MECHANICS_EXTENSION_TESTS Finite rotations, tangent, guides and body force.
names=["finite_rotation_objectivity","nonlinear_consistent_tangent","nonlinear_affine_bar","multipoint_guide","gravity_balance","stable_action","contact_friction_cycle","curved_physical_patch","batched_curved_solid","constrained_modes"];
checks={@objectivity,@tangent,@affine_bar,@guide,@gravity_balance,@stable_action,@contact_cycle,@curved_patch,@batched_curved,@constrained_modes};
metric=zeros(numel(names),1);details=strings(numel(names),1);
for k=1:numel(names),[metric(k),details(k)]=checks{k}();fprintf('%s: true | %s\n',names(k),details(k));end
report=table(names',true(numel(names),1),metric,details,'VariableNames',{'test','passed','metric','details'});
end
function [v,d]=objectivity()
m=fem3d.elevate(fem3d.box_mesh([.1,.02,.01],[1,1,1]));X=m.nodes(m.tetra(1,:),:);a=.7;R=[cos(a),0,sin(a);0,1,0;-sin(a),0,cos(a)];U=X*(R-eye(3))'+[.4,-.2,.3];
[f,~,W]=fem3d.nonlinear_element(X,U,fem3d.isotropic(2e9,.3));v=norm(f);assert(v<1e-6 && W<1e-20);d="40-degree rigid rotation and translation produce zero Green strain and energy";
end
function [v,d]=tangent()
m=fem3d.elevate(fem3d.box_mesh([.1,.02,.01],[1,1,1]));X=m.nodes(m.tetra(1,:),:);rng(41);U=randn(10,3)*1e-4;D=fem3d.isotropic(2e9,.3);[~,K]=fem3d.nonlinear_element(X,U,D);direction=randn(30,1);h=1e-7;
fp=fem3d.nonlinear_element(X,U+h*reshape(direction,3,[])',D);fm=fem3d.nonlinear_element(X,U-h*reshape(direction,3,[])',D);v=norm((fp-fm)/(2*h)-K*direction)/norm(K*direction);
assert(v<1e-7);d="material plus geometric tangent agrees with an independent centered directional derivative";
end
function [v,d]=affine_bar()
m=fem3d.box_mesh([.1,.02,.01],[2,1,1]);mat=struct('D',fem3d.isotropic(2e9,0),'rho',1000);s=fem3d.assemble(m,mat);
X=m.nodes;tip=m.faces(all(reshape(abs(X(m.faces(:),1)-.1)<1e-12,size(m.faces)),2),:);stretch=1.05;P=2e9*.5*(stretch^2-1)*stretch;
F=fem3d.traction(m,tip,[P,0,0]);bc=struct('dofs',fem3d.dofs(find(X(:,1)==0)));r=fem3d.solve_nonlinear(s,F,bc,struct('loadSteps',3));
expected=zeros(size(X));expected(:,1)=(stretch-1)*X(:,1);v=norm(r.displacement-expected,'fro')/norm(expected,'fro');assert(v<1e-6);d="5% finite axial stretch against analytical first Piola load, nu=0";
end
function [v,d]=guide()
m=fem3d.box_mesh([.1,.02,.01],[1,1,1]);s=fem3d.assemble(m,fem3d.material_reference);fixed=fem3d.dofs(find(m.nodes(:,1)==0));tip=find(m.nodes(:,1)==.1);C=sparse(1,3*tip(1),1,1,size(s.K,1));d0=1e-4;
r=fem3d.solve_constrained(s,zeros(size(m.nodes)),struct('dofs',fixed),C,d0);v=abs(r.u(3*tip(1))-d0);assert(v<1e-12 && norm(r.equilibriumResidual(setdiff(1:size(s.K,1),fixed)),inf)<1e-8);
loads=zeros(size(s.K,1),3);loads(3*tip(1),:)=[1,2,3];response=fem3d.constraint_response(s,loads,struct('dofs',fixed),C);free=setdiff(1:size(s.K,1),fixed);oracle=s.K(free,free)\loads(free,:);assert(norm(response.base(free,:)-oracle,'fro')/norm(oracle,'fro')<1e-11);assert(norm(C*response.influence-response.compliance)<1e-12);
d="exact prescribed guide coordinate and reactions; three simultaneous full-size load columns against independent sparse solves";
end
function [v,d]=gravity_balance()
m=fem3d.elevate(fem3d.box_mesh([.1,.02,.01],[2,1,1]));s=fem3d.assemble(m,fem3d.material_reference);F=fem3d.gravity(s,[0,0,-9.81]);v=norm(sum(F,1)-s.totalMass*[0,0,-9.81]);assert(v<1e-12);d="body load integrates the same consistent total mass as dynamics";
end
function [v,d]=stable_action()
m=fem3d.elevate(fem3d.box_mesh([.1,.02,.01],[2,1,1]));s=fem3d.assemble(m,fem3d.material_reference);op=fem3d.action_operator(s);rng(42);U=randn(size(s.K,1),3);v=norm(fem3d.action(op,U)-s.K*U,'fro')/norm(s.K*U,'fro');assert(v<1e-12);d="elementwise strain/stress operator matches assembled stiffness on independent random vectors";
end

function [v,d]=contact_cycle()
s=struct('K',speye(3)*1000);bc=struct('dofs',2);normal=sparse([0,0,1]);tangent=sparse([1,0,0]);previous=struct('tangentialPosition',0);work=0;
for fx=[0:1:6,5:-1:-6,-5:1:0]
 r=fem3d.solve_friction(s,[fx;0;-10],bc,normal,0,tangent,.2,sparse(0,3),zeros(0,1),previous);
 assert(abs(normal*r.u)<1e-12 && abs(r.normalReaction-10)<1e-10 && abs(r.frictionForce)<=2+1e-8);work=work+r.dissipatedWork;previous=r;
end
v=work;assert(work>0);d="unilateral gap, nonnegative reaction, Coulomb stick/slip and dissipative cyclic history";
end
function [v,d]=curved_patch()
A=.2;B=.3;xy=[0,0;1,0;0,1;.5,0;.5,.5;0,.5];m=struct('nodes',[xy,A*xy(:,1).^2+B*xy(:,2).^2],'faces',1:6,'cadFaces',1:6,'cadFaceId',1);
c=struct('geometry',struct('chord_m',1,'span_m',1),'regions',struct('upperTailCadFaces',1),'loads',struct('x_over_c',[.1,.3],'eta',[.6,.7],'integrationDivisions',[4,3]));
p=cad.pressure_patch(m,c);a=.1;b=.3;cc=.1;dd=.2;mx=@(i)(b^(i+1)-a^(i+1))/(i+1);my=@(i)(dd^(i+1)-cc^(i+1))/(i+1);
force=[2*A*mx(1)*my(0),2*B*mx(0)*my(1),-mx(0)*my(0)];
moment=[-mx(0)*my(1)-2*B*(A*mx(2)*my(1)+B*mx(0)*my(3)),mx(1)*my(0)+2*A*(A*mx(3)*my(0)+B*mx(1)*my(2)),2*(B-A)*mx(1)*my(1)];
v=max(norm(sum(p.unitForce,1)-force),norm(sum(cross(m.nodes,p.unitForce,2),1)-moment));assert(v<1e-11 && p.momentDefect_Nm_per_Pa<1e-12);d="fixed physical rectangle on an exact paraboloid Tri6: analytical force, moment and inverse geometry";
end

function [v,d]=batched_curved()
m=fem3d.elevate(fem3d.box_mesh([.1,.02,.01],[1,1,1]));m.nodes(end,3)=m.nodes(end,3)+1e-5;mat=fem3d.material_reference;a=fem3d.assemble(m,mat);m.curvedBatch=true;b=fem3d.assemble(m,mat);rng(55);P=randn(size(a.K,1),3);[K,M]=fem3d.projected_operators(m,mat,P,4);
v=max([norm(a.K-b.K,'fro')/norm(a.K,'fro'),norm(a.M-b.M,'fro')/norm(a.M,'fro'),norm(K-P'*a.K*P,'fro')/norm(K,'fro'),norm(M-P'*a.M*P,'fro')/norm(M,'fro')]);assert(v<1e-11);d="batched curved stiffness/mass and Galerkin operators against independent per-element isoparametric assembly";
end

function [v,d]=constrained_modes()
m=fem3d.box_mesh([.1,.02,.01],[2,1,1]);s=fem3d.assemble(m,fem3d.material_reference);bc=struct('dofs',fem3d.dofs(find(m.nodes(:,1)==0)));tip=find(m.nodes(:,1)==.1);C=sparse(1,3*tip(1),1,1,size(s.K,1));r=fem3d.modal_constrained(s,bc,C,3);free=setdiff(1:size(s.K,1),bc.dofs);Z=null(full(C(:,free)));[~,D]=eig(Z'*s.K(free,free)*Z,Z'*s.M(free,free)*Z);exact=sort(sqrt(diag(D)))/(2*pi);v=max(abs(r.frequencyHz./exact(1:3)-1));assert(v<1e-8 && max(r.constraintResidual)<1e-10 && max(r.eigenEnergyResidual)<1e-7);d="projected stiffness inverse eigenmodes against independent explicit nullspace generalized eigenproblem";
end
