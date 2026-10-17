function patch = curved_patch(mesh,config,use)
%CURVED_PATCH Integrate a fixed physical rectangle on curved graph skins.
% Composite positive Gauss quadrature in physical x/y; inverse FE mapping is
% evaluated at each point. Projected area is exact and independent of the mesh.
% Refining integrationDivisions controls the curved/element-interface integral.
box=[config.loads.x_over_c(:)'*config.geometry.chord_m;(config.loads.eta(:)'-.5)*config.geometry.span_m];
divisions=[32,16];if isfield(config.loads,'integrationDivisions'),divisions=config.loads.integrationDivisions(:)';end
[x,w]=fem3d.gauss01(2);nq=4*prod(divisions);I=zeros(nq,6);J=I;values=I;points=zeros(nq,3);normals=points;weights=zeros(nq,1);owners=zeros(nq,1);row=0;
delta=diff(box,1,2)./divisions';
for ix=1:divisions(1),for iy=1:divisions(2),for a=1:2,for b=1:2
    xy=box(:,1)'+([ix-1,iy-1]+[x(a),x(b)]).*delta';hit=cad.ray(mesh,xy,use);
    assert(~isempty(hit.z),'cad:Patch','Curved physical patch has a coverage gap.');
    f=hit.face(end);face=mesh.faces(f,:);r=hit.parametric(end,:);[N,G]=fem3d.tri_shape(r,6);X=mesh.nodes(face,:);
    tangent=X'*G;normal=cross(tangent(:,1),tangent(:,2))';normal=normal/norm(normal);
    assert(normal(3)>1e-5,'cad:Patch','Physical rectangle needs a regular upward graph skin.');
    row=row+1;I(row,:)=row;J(row,:)=face;values(row,:)=N;points(row,:)=N*X;normals(row,:)=normal;
    weights(row)=prod(delta)*w(a)*w(b)/normal(3);owners(row)=f;
end,end,end,end
H=sparse(I(:),J(:),values(:),nq,size(mesh.nodes,1));P=-weights.*normals;F=H'*P;
forceDefect=norm(sum(F,1)-sum(P,1));momentDefect=norm(sum(cross(mesh.nodes,F,2),1)-sum(cross(points,P,2),1));
projected=sum(weights.*normals(:,3));expected=prod(diff(box,1,2));
assert(abs(projected/expected-1)<1e-10,'cad:Patch','Projected physical quadrature area mismatch.');
patch=struct('H',H,'points',points,'weights',weights,'normals',normals,'faces',owners,'unitForce',F, ...
    'area_m2',sum(weights),'projectedArea_m2',projected,'expectedProjectedArea_m2',expected, ...
    'forceDefect_N_per_Pa',forceDefect,'momentDefect_Nm_per_Pa',momentDefect,'bounds',box, ...
    'pressureSign','positive inward','integrationDivisions',divisions, ...
    'method','composite physical-rectangle Gauss quadrature with exact inverse Tri6 mapping');
% Independent exposure checks at points not tied to the quadrature grid.
for xx=linspace(box(1,1),box(1,2),5),for yy=linspace(box(2,1),box(2,2),3)
    allHit=cad.ray(mesh,[xx,yy]);part=cad.ray(mesh,[xx,yy],use);
    assert(~isempty(part.z) && abs(part.z(end)-allHit.z(end))<1e-9,'cad:Region','Curved patch is not exposed to the exterior.');
end,end
end
