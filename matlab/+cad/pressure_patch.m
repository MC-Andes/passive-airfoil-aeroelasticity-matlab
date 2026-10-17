function patch = pressure_patch(mesh,config)
%PRESSURE_PATCH Clip explicit exterior CAD skins by a fixed physical rectangle.
% Straight Tri3/Tri6 only. Quadrature on clipped subtriangles maps back to
% parent FE shape functions, so load footprint is independent of mesh size.
labels=cad.face_labels(mesh); use=find(ismember(labels,config.regions.upperTailCadFaces));
assert(~isempty(use),'cad:Region','Missing explicitly identified upper tail surfaces.');
if size(mesh.faces,2)==6
    F=mesh.faces(use,:);X=mesh.nodes;curved=false;
    for k=1:3
        a=k;b=mod(k,3)+1;
        if max(abs(X(F(:,k+3),:)-(X(F(:,a),:)+X(F(:,b),:))/2),[],'all')>1e-10,curved=true;break;end
    end
    if curved,patch=cad.curved_patch(mesh,config,use);return;end
end
box=[config.loads.x_over_c(:)'*config.geometry.chord_m; ...
    (config.loads.eta(:)'-.5)*config.geometry.span_m];
rows={}; cols={}; vals={}; points=[]; normals=[]; weights=[]; owners=[]; row=0;
for f=use'
    face=mesh.faces(f,:); xe=mesh.nodes(face,:);
    if size(face,2)==6
        assert(max(abs(xe(4:6,:)-(xe([1,2,3],:)+xe([2,3,1],:))/2),[],'all')<1e-10, ...
            'cad:Patch','Curved-face clipping needs parametric integration; use full-face quadrature.');
    end
    % Rows carry physical coordinates and barycentric coordinates in the parent.
    poly=[xe(1:3,:),eye(3)];
    for axis=1:2
        poly=clip(poly,axis,box(axis,1),1); poly=clip(poly,axis,box(axis,2),-1);
    end
    if size(poly,1)<3, continue; end
    n=cross(xe(2,:)-xe(1,:),xe(3,:)-xe(1,:)); n=n/norm(n);
    assert(n(3)>0,'cad:Region','Selected upper skin has inward/downward orientation.');
    for j=2:size(poly,1)-1
        tri=poly([1,j,j+1],:); area=norm(cross(tri(2,1:3)-tri(1,1:3),tri(3,1:3)-tri(1,1:3)))/2;
        if area<1e-20, continue; end
        % Degree-2 rule integrates Tri6 nodal forces and moment on a flat subtriangle.
        bary=[2/3,1/6,1/6;1/6,2/3,1/6;1/6,1/6,2/3]*tri(:,4:6);
        for k=1:3
            row=row+1;
            if numel(face)==3, sh=bary(k,:); else, sh=fem3d.tri_shape(bary(k,2:3),6); end
            rows{row}=repmat(row,numel(face),1); cols{row}=face'; vals{row}=sh'; %#ok<AGROW>
            points(row,:)=sh*xe; normals(row,:)=n; weights(row,1)=area/3; owners(row,1)=f; %#ok<AGROW>
        end
    end
end
assert(row>0,'cad:Patch','Physical patch has no surface coverage.');
H=sparse(vertcat(rows{:}),vertcat(cols{:}),vertcat(vals{:}),row,size(mesh.nodes,1));
projected=sum(weights.*normals(:,3)); expected=prod(diff(box,1,2));
assert(abs(projected/expected-1)<1e-6,'cad:Patch','Patch has gaps or overlapping skins; inspect region labels.');
% Verify exposed upper surface at a fixed grid independent of the structural mesh.
xx=linspace(box(1,1),box(1,2),7); yy=linspace(box(2,1),box(2,2),5);
for x=xx, for y=yy
    allHit=cad.ray(mesh,[x,y]); selected=cad.ray(mesh,[x,y],use);
    assert(~isempty(selected.z) && ~isempty(allHit.z) && abs(selected.z(end)-allHit.z(end))<1e-9, ...
        'cad:Region','Pressure patch does not lie on the exposed upper envelope.');
end, end
unitForce=H'*(-weights.*normals); P=-weights.*normals;
forceDefect=norm(sum(unitForce,1)-sum(P,1));
momentDefect=norm(sum(cross(mesh.nodes,unitForce,2),1)-sum(cross(points,P,2),1));
patch=struct('H',H,'points',points,'weights',weights,'normals',normals,'faces',owners, ...
    'unitForce',unitForce,'area_m2',sum(weights),'projectedArea_m2',projected, ...
    'expectedProjectedArea_m2',expected,'forceDefect_N_per_Pa',forceDefect, ...
    'momentDefect_Nm_per_Pa',momentDefect,'bounds',box,'pressureSign','positive inward');
end

function out=clip(poly,axis,bound,sgn)
out=zeros(0,6); if isempty(poly), return; end
for k=1:size(poly,1)
    a=poly(k,:); b=poly(mod(k,size(poly,1))+1,:); da=sgn*(a(axis)-bound); db=sgn*(b(axis)-bound);
    if da>=0, out(end+1,:)=a; end %#ok<AGROW>
    if (da>=0)~=(db>=0), out(end+1,:)=a+da/(da-db)*(b-a); end %#ok<AGROW>
end
end
