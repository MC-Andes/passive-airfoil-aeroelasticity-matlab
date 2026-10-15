function results = run_fem_tests()
%RUN_FEM_TESTS Independent analytical/conservation checks of own FEM core.
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(fullfile(root,'matlab'));
names=["tetra_rigid_body","consistent_mass","orientation_invariance", ...
    "affine_patch_test","pressure_balance","cantilever_refinement", ...
    "modal_residual","newmark_free_mode","conservative_transfer", ...
    "material_orthotropic_limit","unsupported_mesh_rejected","cad_volume"];
checks={@tetra_rigid_body,@consistent_mass,@orientation_invariance,@affine_patch_test, ...
    @pressure_balance,@cantilever_refinement,@modal_residual,@newmark_free_mode, ...
    @conservative_transfer,@material_orthotropic_limit,@unsupported_mesh_rejected,@cad_volume};
passed=false(numel(names),1); metric=NaN(numel(names),1); details=strings(numel(names),1);
for k=1:numel(names)
    try
        [metric(k),details(k)]=checks{k}(); passed(k)=true;
    catch err
        details(k)=string(err.identifier)+": "+string(err.message);
    end
    fprintf('%s: %s | %s\n',names(k),string(passed(k)),details(k));
end
results=table(names',passed,metric,details,'VariableNames',{'test','passed','metric','details'});
output=fullfile(root,'results'); if ~exist(output,'dir'), mkdir(output); end
writetable(results,fullfile(output,'fem_test_results.csv'));
assert(all(passed),'fem3d:Tests','One or more FEM verification tests failed.');
end

function [metric,detail]=tetra_rigid_body()
X=[0,0,0;2,0,0;.1,1,0;.2,.1,.7]; D=fem3d.isotropic(3e9,.3);
[K,~,B]=fem3d.tet4(X,D,1200); U=[.2,-.4,.1]+cross(repmat([.1,.2,-.3],4,1),X,2);
u=reshape(U',[],1); metric=norm(B*u);
ev=sort(eig(K)); count=sum(abs(ev)<1e-9*norm(K,2));
assert(metric<1e-13 && count==6 && ev(7)>0); detail="six rigid modes; zero strain under rigid motion";
end

function [metric,detail]=consistent_mass()
X=[0,0,0;1,0,0;0,1,0;0,0,1]; rho=1200;
[~,M,~,vol]=fem3d.tet4(X,eye(6),rho); v=repmat([.3;-.2;.7],4,1);
expected=rho*vol*sum([.3,-.2,.7].^2); metric=abs(v'*M*v-expected)/expected;
assert(metric<1e-13 && min(eig(M))>0); detail="translational kinetic energy equals rho*V*|velocity|^2";
end

function [metric,detail]=orientation_invariance()
mesh=fem3d.box_mesh([.2,.1,.05],[2,2,2]); mat=fem3d.material_reference;
a=fem3d.assemble(mesh,mat); mesh.tetra(:,[2,3])=mesh.tetra(:,[3,2]);
b=fem3d.assemble(mesh,mat); metric=norm(a.K-b.K,'fro')/norm(a.K,'fro');
assert(metric<1e-13 && norm(a.M-b.M,'fro')<1e-12); detail="connectivity orientation does not alter assembled physics";
end

function [metric,detail]=affine_patch_test()
mesh=fem3d.box_mesh([.2,.12,.08],[3,3,3]); mat=fem3d.material_reference;
system=fem3d.assemble(mesh,mat); A=[1,.2,-.3;.1,-.4,.2;.05,-.1,.3]*1e-4;
U=mesh.nodes*A'+[2e-5,-1e-5,3e-5]; expected=reshape(U',[],1);
epsv=[A(1,1);A(2,2);A(3,3);A(1,2)+A(2,1);A(2,3)+A(3,2);A(1,3)+A(3,1)];
s=mat.D*epsv; S=[s(1),s(4),s(6);s(4),s(2),s(5);s(6),s(5),s(3)];
F=mesh.faces; X=mesh.nodes;
n=cross(X(F(:,2),:)-X(F(:,1),:),X(F(:,3),:)-X(F(:,1),:),2);
n=n./vecnorm(n,2,2);
loads=fem3d.traction(mesh,F,n*S');
fixed=fem3d.dofs(find(X(:,1)==0)); bc=struct('dofs',fixed,'values',expected(fixed));
r=fem3d.solve_static(system,loads,bc); metric=norm(r.u-expected,inf)/norm(expected,inf);
assert(metric<1e-10 && max(abs(r.stress.strain-epsv'),[],'all')<1e-12);
assert(norm(r.forceBalance)<1e-8 && norm(r.momentBalance)<1e-8);
detail="exact affine solution under analytically derived boundary tractions and nonzero BC";
end

function [metric,detail]=pressure_balance()
mesh=fem3d.box_mesh([.2,.12,.08],[3,2,2]); F=fem3d.traction(mesh,mesh.faces,1000);
metric=norm(sum(F,1)); moment=norm(sum(cross(mesh.nodes,F,2),1));
assert(metric<1e-12 && moment<1e-12); detail="uniform pressure on closed surface has zero force and moment";
end

function [metric,detail]=cantilever_refinement()
E=2617e6; nu=.3; mat=struct('D',fem3d.isotropic(E,nu),'rho',1160);
L=.2; width=.04; h=.02; analytical=L^3/(3*E*(width*h^3/12));
divisions=[8,2,2;16,4,4;32,8,8]; errors=zeros(3,1); deflection=errors;
for k=1:3
    mesh=fem3d.box_mesh([L,width,h],divisions(k,:)); system=fem3d.assemble(mesh,mat);
    X=mesh.nodes;
    endfaces=mesh.faces(all(reshape(abs(X(mesh.faces(:),1)-L)<1e-12,size(mesh.faces)),2),:);
    F=fem3d.traction(mesh,endfaces,[0,0,-1/(width*h)]);
    bc=struct('dofs',fem3d.dofs(find(X(:,1)==0)));
    result=fem3d.solve_static(system,F,bc);
    tip=find(abs(X(:,1)-L)<1e-12);
    deflection(k)=-mean(result.displacement(tip,3)); errors(k)=abs(deflection(k)/analytical-1);
end
metric=errors(end); assert(all(diff(errors)<0) && metric<.15);
detail=sprintf('Euler-Bernoulli reference; relative errors coarse/mid/fine = %.4g / %.4g / %.4g',errors);
end

function [system,bc,modes]=modal_case()
mesh=fem3d.box_mesh([.1,.02,.01],[12,3,3]); mat=fem3d.material_reference;
system=fem3d.assemble(mesh,mat); bc=struct('dofs',fem3d.dofs(find(mesh.nodes(:,1)==0)));
modes=fem3d.modal(system,bc,3);
end

function [metric,detail]=modal_residual()
[system,~,modes]=modal_case(); metric=max(modes.eigenResidual);
assert(metric<1e-8 && norm(modes.Phi'*system.M*modes.Phi-eye(3),'fro')<1e-8);
detail="assembled generalized eigenproblem and mass orthogonality";
end

function [metric,detail]=newmark_free_mode()
[system,bc,modes]=modal_case(); u0=modes.Phi(:,1); u0=u0*1e-4/max(abs(u0));
T=2*pi/modes.omega(1); time=linspace(0,2*T,801);
r=fem3d.newmark(system,bc,time,@(t) zeros(size(u0)),u0);
expected=u0*cos(modes.omega(1)*time);
metric=max(abs(r.U-expected),[],'all')/max(abs(u0));
energyError=max(abs(r.energy/r.energy(1)-1));
assert(metric<5e-4 && energyError<1e-8 && max(r.freeResidual)<1e-6);
detail=sprintf('exact cosine first mode; amplitude error %.3g; energy drift %.3g',metric,energyError);
end

function [metric,detail]=conservative_transfer()
rng(17); X=rand(100,3); P=rand(7,3); force=randn(7,3);
[H,~]=fem3d.transfer(X,P,(1:100)',16);
F=H'*force; A=randn(3); U=X*A'+randn(1,3);
affineError=norm(H*X-P,'fro');
forceError=norm(sum(F,1)-sum(force,1));
momentError=norm(sum(cross(X,F,2),1)-sum(cross(P,force,2),1));
workError=abs(sum(U.*F,'all')-sum((H*U).*force,'all'));
metric=max([affineError,forceError,momentError,workError]);
assert(metric<1e-12); detail="affine kinematics, force, moment and virtual work conserved";
end

function [metric,detail]=material_orthotropic_limit()
E=2.5e9; nu=.32; G=E/(2*(1+nu));
D=fem3d.orthotropic([E,E,E],[nu,nu,nu],[G,G,G]); I=fem3d.isotropic(E,nu);
metric=norm(D-I,'fro')/norm(I,'fro'); assert(metric<1e-13);
detail="orthotropic compliance reduces to isotropic constitutive law";
end

function [metric,detail]=unsupported_mesh_rejected()
mesh=fem3d.box_mesh([.1,.1,.1],[1,1,1]); mat=fem3d.material_reference;
s=fem3d.assemble(mesh,mat); F=zeros(size(mesh.nodes)); F(8,1)=1;
rejected=false;
try
    fem3d.solve_static(s,F,struct('dofs',1));
catch e
    rejected=strcmp(e.identifier,'fem3d:Unconstrained');
end
assert(rejected); metric=1; detail="underconstrained solid rejected instead of reporting false equilibrium";
end

function [metric,detail]=cad_volume()
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
mesh=fem3d.read_mesh(fullfile(root,'data','fem_reference'));
T=mesh.tetra; X=mesh.nodes; volume=0;
for k=1:size(T,1)
    a=X(T(k,:),:); volume=volume+abs(det((a(2:4,:)-a(1,:))'))/6;
end
integrationError=abs(volume-mesh.info.total_volume_m3)/volume;
audit=jsondecode(fileread(fullfile(root,'data','reference','cad_volume_audit.json')));
metric=abs(volume-audit.cad_volume_m3)/audit.cad_volume_m3;
assert(integrationError<1e-10 && metric<.01 && audit.brep_valid && ...
    strcmp(audit.source_sha256,mesh.info.source_sha256) && size(mesh.nodes,1)==mesh.info.node_count && ...
    size(mesh.faces,1)==mesh.info.face_count);
detail=sprintf('volume differs %.4g from independently integrated CAD; valid BREP and exterior faces',metric);
end
