function results = run_cad_tests()
%RUN_CAD_TESTS Quadratic solids and CAD operators; independent analytical checks.
names=["tet10_rigid_mass","tet10_affine_patch","tet10_permutation", ...
    "tet10_curved_quadrature","tet10_fold_rejected","tri6_curved_load", ...
    "tet10_bending_modes","axis_reflection","physical_patch_and_sensors","tet10_transient","provenance_hash"];
checks={@rigid_mass,@affine_patch,@permutation,@curved_quadrature,@fold_rejected, ...
    @curved_load,@bending_modes,@axis_reflection,@patch_sensors,@tet10_transient,@provenance_hash};
passed=false(numel(names),1); metric=NaN(size(passed)); details=strings(size(passed));
for k=1:numel(names)
    try, [metric(k),details(k)]=checks{k}(); passed(k)=true;
    catch err, details(k)=string(err.identifier)+": "+string(err.message); end
    fprintf('%s: %s | %s\n',names(k),string(passed(k)),details(k));
end
results=table(names',passed,metric,details,'VariableNames',{'test','passed','metric','details'});
assert(all(passed),'cad:Tests','CAD/Tet10 verification failed.');
end

function [metric,detail]=rigid_mass()
m=fem3d.elevate(fem3d.box_mesh([1,1,1],[1,1,1])); X=m.nodes(m.tetra(1,:),:);
[K,M,B,V]=fem3d.tet10(X,fem3d.isotropic(3e9,.3),1200);
U=[.1,-.2,.4]+cross(repmat([.2,.1,-.3],10,1),X,2); u=reshape(U',[],1);
ev=sort(eig(K)); assert(sum(abs(ev)<1e-9*norm(K,2))==6 && ev(7)>0);
v=repmat([.3;-.1;.6],10,1); metric=abs(v'*M*v/(1200*V*sum([.3,-.1,.6].^2))-1);
assert(metric<1e-12 && norm(B*u)<1e-14 && norm(K*u)<1e-5 && min(eig(M))>0);
detail="six rigid motions, exact translational kinetic energy, positive consistent mass";
end

function [metric,detail]=affine_patch()
m=fem3d.elevate(fem3d.box_mesh([.2,.12,.08],[2,2,2])); mat=fem3d.material_reference;
s=fem3d.assemble(m,mat); A=[1,.2,-.3;.1,-.4,.2;.05,-.1,.3]*1e-4;
U=m.nodes*A'+[2e-5,-1e-5,3e-5]; u=reshape(U',[],1);
e=[A(1,1);A(2,2);A(3,3);A(1,2)+A(2,1);A(2,3)+A(3,2);A(1,3)+A(3,1)];
sig=mat.D*e; S=[sig(1),sig(4),sig(6);sig(4),sig(2),sig(5);sig(6),sig(5),sig(3)];
F=fem3d.traction(m,m.faces,@(p,n)n*S');
fixed=fem3d.dofs(find(m.nodes(:,1)==0)); r=fem3d.solve_static(s,F,struct('dofs',fixed,'values',u(fixed)));
metric=norm(r.u-u,inf)/norm(u,inf);
assert(metric<1e-10 && max(abs(r.stress.strain-e'),[],'all')<1e-12);
q=fem3d.surface_quadrature(m,m.faces); P=q.weights.*(q.normals*S');
assert(norm(q.H*m.nodes-q.points,'fro')<1e-13 && ...
    abs(sum((q.H*U).*P,'all')-sum(U.*(q.H'*P),'all'))<1e-12);
detail="exact affine solid with nonzero supports and Tri6 analytical boundary tractions; adjoint work";
end

function [metric,detail]=permutation()
m=fem3d.elevate(fem3d.box_mesh([.2,.1,.05],[2,1,1])); mat=fem3d.material_reference;
a=fem3d.assemble(m,mat);
% Swap vertices 2/3 and the associated edges, including a reflection.
m.tetra=m.tetra(:,[1,3,2,4,7,6,5,8,10,9]);
[m.faces,m.faceOwners]=fem3d.boundary(m.tetra,m.nodes); b=fem3d.assemble(m,mat);
metric=norm(a.K-b.K,'fro')/norm(a.K,'fro');
assert(metric<1e-12 && norm(a.M-b.M,'fro')/norm(a.M,'fro')<1e-12);
assert(norm(sum(fem3d.traction(m,m.faces,100),1))<1e-12);
detail="vertex/edge permutation, orientation and closed Tri6 pressure balance";
end

function [metric,detail]=curved_quadrature()
m=fem3d.elevate(fem3d.box_mesh([1,1,1],[1,1,1])); X=m.nodes(m.tetra(1,:),:);
X(:,3)=X(:,3)+.08*X(:,1).*X(:,2); D=fem3d.isotropic(3e9,.3);
[K4,M4,~,V4]=fem3d.tet10(X,D,1200,4); [K6,M6,~,V6]=fem3d.tet10(X,D,1200,6);
[K8,~,~,V8]=fem3d.tet10(X,D,1200,8);
metric=norm(K6-K8,'fro')/norm(K8,'fro');
assert(metric<1e-9 && norm(K4-K8,'fro')/norm(K8,'fro')<1e-6 && ...
    norm(M4-M6,'fro')/norm(M6,'fro')<1e-10 && abs(V4-V6)<1e-12 && abs(V6-V8)<1e-12);
U=X*[.1,.02,.03;0,.2,0;0,.01,.3]'; [~,~,B]=fem3d.tet10(X,D,1200);
assert(norm(B*reshape(U',[],1)-[.1;.2;.3;.02;.01;.03])<1e-13);
detail="curved isoparametric mapping; independent 4/6/8-point-per-axis integration and affine reproduction";
end

function [metric,detail]=fold_rejected()
m=fem3d.elevate(fem3d.box_mesh([1,1,1],[1,1,1])); X=m.nodes(m.tetra(1,:),:); X(5,:)=X(5,:)+[0,4,3];
rejected=false; try, fem3d.tet10(X,eye(6),1); catch e, rejected=ismember(e.identifier,{'fem3d:Jacobian','fem3d:Degenerate'}); end
assert(rejected); metric=1; detail="folded quadratic geometry rejected at extra control/quadrature points";
end

function [metric,detail]=curved_load()
xy=[0,0;1,0;0,1;.5,0;.5,.5;0,.5]; m.nodes=[xy,.2*xy(:,1).^2+.1*xy(:,2).^2]; f=1:6;
q=fem3d.surface_quadrature(m,f,5); F=fem3d.traction(m,f,7);
% Graph normal dA = [-.4*x,-.2*y,1] dxdy; triangle area=1/2, first moments=1/6.
expected=7*[.4/6,.2/6,-.5]; metric=norm(sum(F,1)-expected);
P=-7*q.weights.*q.normals;
assert(metric<1e-12 && norm(sum(cross(m.nodes,F,2),1)-sum(cross(q.points,P,2),1))<1e-12);
detail="curved Tri6 pressure integrates exact graph resultant and independent higher-order moment";
end

function [metric,detail]=bending_modes()
L=.2; b=.02; h=.004; E=2.6e9; rho=1160; mat=struct('D',fem3d.isotropic(E,0),'rho',rho);
I=b*h^3/12; exact=L^3/(3*E*I); fexact=1.875104068711961^2/(2*pi*L^2)*sqrt(E*I/(rho*b*h));
err=zeros(3,1); ferr=err;
for k=1:3
    m=fem3d.elevate(fem3d.box_mesh([L,b,h],[4*2^(k-1),2,1])); s=fem3d.assemble(m,mat);
    mask=all(reshape(abs(m.nodes(m.faces(:),1)-L)<1e-12,size(m.faces)),2); faces=m.faces(mask,:);
    F=fem3d.traction(m,faces,[0,0,-1/(b*h)]); bc=struct('dofs',fem3d.dofs(find(m.nodes(:,1)==0)));
    r=fem3d.solve_static(s,F,bc); def=-mean(r.displacement(abs(m.nodes(:,1)-L)<1e-12,3));
    modes=fem3d.modal(s,bc,2); err(k)=abs(def/exact-1); ferr(k)=abs(modes.frequencyHz(1)/fexact-1);
    assert(max(modes.eigenResidual)<1e-8);
end
metric=err(end); assert(all(diff(err)<0) && err(end)<.02 && ferr(end)<.01);
detail=sprintf('slender beam analytical bending errors %.3g %.3g %.3g; frequency errors %.3g %.3g %.3g',err,ferr);
end

function [metric,detail]=axis_reflection()
t=cad.transform([1,3,2],[.1,.2,.3]); x=[.2;.4;.6]; F=[1;3;2]; R=t.R;
metric=norm(cross(R*x,R*F)-t.momentMap*cross(x,F));
assert(t.determinant==-1 && metric<1e-14);
detail="axial moments under improper CAD axis permutation verified by r cross F";
end

function [metric,detail]=patch_sensors()
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); c=jsondecode(fileread(fullfile(root,'cases','cad_original','case.json')));
m=fem3d.read_mesh(fullfile(root,'data','fem_reference'),c.geometry.axisOrder,c.geometry.origin_m);
p=cad.pressure_patch(m,c); s=cad.sensors(m,c);
U=m.nodes*[.1,.02,0;.03,.2,0;.1,-.1,.3]'+[.2,.1,-.2]; F=p.unitForce;
metric=max([p.forceDefect_N_per_Pa,p.momentDefect_Nm_per_Pa, ...
    abs(sum(U.*F,'all')-sum((p.H*U).*(-p.weights.*p.normals),'all'))]);
assert(metric<1e-12 && norm(s.H*m.nodes-s.points,'fro')<1e-12);
assert(norm(sum(F,1).*[0,0,1]-[0,0,-p.expectedProjectedArea_m2])<1e-12);
detail="source CAD exposed fixed patch, exact projected area and force/moment/work; nine material sensors";
end

function [metric,detail]=tet10_transient()
m=fem3d.elevate(fem3d.box_mesh([.1,.02,.01],[4,1,1])); s=fem3d.assemble(m,fem3d.material_reference);
bc=struct('dofs',fem3d.dofs(find(m.nodes(:,1)==0))); mo=fem3d.modal(s,bc,2);
u0=mo.Phi(:,1)*1e-5/max(abs(mo.Phi(:,1))); time=linspace(0,2/mo.frequencyHz(1),401);
r=fem3d.newmark(s,bc,time,@(t)zeros(size(u0)),u0);
metric=max(abs(r.U-u0*cos(mo.omega(1)*time)),[],'all')/max(abs(u0));
assert(metric<.001 && max(abs(r.energy/r.energy(1)-1))<1e-8);
detail="Tet10 free vibration against exact modal cosine; Newmark energy conservation";
end

function [metric,detail]=provenance_hash()
p=[tempname,"_a' b.txt"]; p=char(join(p,''));
c=onCleanup(@()delete(p)); fid=fopen(p,'wb'); fwrite(fid,uint8('abc'),'uint8'); fclose(fid);
value=cad.sha256(p); assert(strcmp(value,'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad'));
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); config=cad.load_case(root);
assert(config.flow.Re>1e5 && numel(config.caseSha256)==64); metric=0;
detail="known SHA256 vector, quoted path with spaces/apostrophe, versioned source-audited case";
end
