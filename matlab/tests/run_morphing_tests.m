function report = run_morphing_tests()
%RUN_MORPHING_TESTS Updated study references, physical supports and actuation.
names=["morphing_flow_reference","shaft_moment_convention","displacement_control_modulus","bonded_socket_and_servo","nonlinear_batch_oracle","zero_force_finite_actuation"];
checks={@flow_reference,@coefficient_moment,@modulus,@physical_bc,@batch_oracle,@finite_command};
passed=false(numel(names),1);error=passed+NaN;description=strings(size(names));
for k=1:numel(names),[error(k),description(k)]=checks{k}();passed(k)=true;fprintf('%s: true | %s\n',names(k),description(k));end
report=table(names(:),passed,error,description(:),'VariableNames',{'test','passed','error','description'});
end
function c=config()
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));c=morphing.load_case(root);
end
function [v,d]=flow_reference()
c=config();q=[];Re=[];for U=c.flow.speeds_m_s(:)',r=morphing.aerodynamic_coefficients(zeros(0,3),zeros(0,3),c,U,0);q(end+1)=r.q_Pa;Re(end+1)=r.Re;end
v=max(abs(q-[99,275,539]));assert(v<1e-12 && max(abs(Re./[1.84e5,3.07e5,4.30e5]-1))<.003 && abs(2*c.geometry.span_m/c.geometry.chord_m-2.88)<1e-12);
d="Bogota rho=.88, nominal c=.25 and half-wing S=.09; q and Re independently checked";
end
function [v,d]=coefficient_moment()
c=config();U=25;a=30;points=[.2,.1,.02;.1,.3,-.01];F=[3,1,4;-1,2,-1];r=morphing.aerodynamic_coefficients(points,F,c,U,a);
% Independent component expression about actual offset shaft, not LE/zero-height.
my=sum((points(:,3)-.0043).*F(:,1)-(points(:,1)-.0625).*F(:,3));q=.5*.88*25^2;
v=max(abs([r.moment_Nm(2)-my,r.CD-(2*cosd(a)+3*sind(a))/(q*.09),r.CL-(-2*sind(a)+3*cosd(a))/(q*.09),r.Cm-my/(q*.09*.25)]));assert(v<1e-13);
d="actual force arms, wind axes and pitching moment about c/4 plus 4.3 mm shaft offset";
end
function [v,d]=modulus()
m=fem3d.box_mesh([.1,.02,.01],[2,1,1]);E=1.6e9;mat=struct('E',E,'nu',0,'rho',1270,'D',fem3d.isotropic(E,0));s=fem3d.assemble(m,mat);tip=find(abs(m.nodes(:,1)-.1)<1e-12);base=find(m.nodes(:,1)<1e-12);delta=1e-5;bc=struct('dofs',[fem3d.dofs(base);fem3d.dofs(tip,1)],'values',[zeros(3*numel(base),1);repmat(delta,numel(tip),1)]);
r=fem3d.solve_static(s,zeros(size(s.K,1),1),bc);s2=s;s2.K=2*s.K;r2=fem3d.solve_static(s2,zeros(size(s.K,1),1),bc);
v=norm(r.u-r2.u,inf);assert(v<1e-15 && abs(sum(r.reaction(tip,1))/(E*.02*.01*delta/.1)-1)<1e-10 && norm(r.displacement(:,1)-delta*m.nodes(:,1)/.1,inf)<1e-14);
d="displacement-driven homogeneous bar: exact geometry independent of uniform E, reaction proportional to E; no spurious load-control limit";
end
function [v,d]=physical_bc()
c=config();h=c.audit.shaft_socket.side_m/2;
a=fem3d.box_mesh([2*h,.36,2*h],[1,4,1]);a.nodes=a.nodes+[.0625-h,0,.0043-h];
x=c.audit.slider_drive.native_plane_x_m-c.geometry.origin_m(1);b=fem3d.box_mesh([.01,.36,.003],[1,4,1]);b.nodes=b.nodes+[x,0,-.006];n=size(a.nodes,1);m=struct('nodes',[a.nodes;b.nodes],'tetra',[a.tetra;b.tetra+n]);[m.faces,~]=fem3d.boundary(m.tetra,m.nodes);m.cadFaces=m.faces;m.cadFaceId=ones(size(m.faces,1),1)*999;m.info=struct('source_sha256',c.geometry.sourceSha256);
for k=1:size(m.faces,1)
 X=m.nodes(m.faces(k,:),:);
 if all(abs(X(:,1)-(.0625-h))<1e-10),m.cadFaceId(k)=74;
 elseif all(abs(X(:,1)-(.0625+h))<1e-10),m.cadFaceId(k)=72;
 elseif all(abs(X(:,3)-(.0043-h))<1e-10),m.cadFaceId(k)=73;
 elseif all(abs(X(:,3)-(.0043+h))<1e-10),m.cadFaceId(k)=71;
 elseif all(abs(X(:,1)-x)<1e-10),m.cadFaceId(k)=8;
 end
end
[bc,info]=morphing.boundary_conditions(m,c,-1e-5);u=zeros(3*size(m.nodes,1),1);u(bc.dofs)=bc.values;v=max(abs(u(fem3d.dofs(info.driveNodes,1))+1e-5));
assert(v<1e-15 && all(m.nodes(info.mountNodes,2)<=.324+1e-10) && ~any(ismember(fem3d.dofs(find(abs(m.nodes(:,2)-.36)<1e-10 & (1:size(m.nodes,1))'<=n)),bc.dofs)));
bad=m;bad.info.source_sha256='historical';failed=false;try,morphing.boundary_conditions(bad,c,0);catch err,failed=strcmp(err.identifier,'morphing:Hash');end;assert(failed);
d="spatially distinct bonded socket/servo, insertion truncation, free shaft-tip region and rejection of historical face labels";
end
function [v,d]=batch_oracle()
rng(312);X=repmat([0,0,0;.12,0,0;0,.09,0;0,0,.07],1,1,5);U=.003*randn(size(X));D=fem3d.isotropic(1.6e9,.35);[f,K,W,J]=fem3d.nonlinear_tet4_batch(X,U,D);v=0;referenceJ=zeros(5,1);
for k=1:5,[ff,KK,ww,jj]=fem3d.nonlinear_element(X(:,:,k),U(:,:,k),D);v=max([v,norm(f(:,k)-ff)/norm(ff),norm(K(:,:,k)-KK,'fro')/norm(KK,'fro'),abs(W(k)/ww-1)]);referenceJ(k)=jj;end
v=max(v,abs(J-min(referenceJ)));
assert(v<1e-12);d="batched finite-strain force, tangent and energy against independent per-element SVK implementation on random states";
end
function [v,d]=finite_command()
m=fem3d.box_mesh([.1,.02,.01],[2,1,1]);m.nonlinearBatch=true;mat=struct('D',fem3d.isotropic(1.6e9,0),'rho',1270);s=fem3d.assemble(m,mat);tip=find(m.nodes(:,1)>.099);base=find(m.nodes(:,1)<1e-12);delta=.005;bc=struct('dofs',[fem3d.dofs(base);fem3d.dofs(tip,1)],'values',[zeros(3*numel(base),1);repmat(delta,numel(tip),1)]);
r=fem3d.solve_nonlinear(s,zeros(size(s.K,1),1),bc,struct('loadSteps',4));v=norm(r.displacement(:,1)-delta*m.nodes(:,1)/.1,inf);lambda=1.05;P=1.6e9*(lambda^2-1)*lambda/2;force=reshape(r.state.force,3,[])';assert(v<1e-11 && abs(sum(force(tip,1))/(P*.02*.01)-1)<1e-7);
d="zero external-force, prescribed 5% finite bar stretch: exact first-Piola reaction and convergence normalized to driven reactions";
end
