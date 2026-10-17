function transfer = section_traction_transfer(mesh,config,rows,Re,alpha)
%SECTION_TRACTION_TRANSFER Extrude computed 2D skin tractions to original skins.
% rows come from nsflow.surface_traction. This is explicitly a span-uniform
% 2D loading hypothesis; no finite-span pressure claim is implied.
chord=config.geometry.chord_m;span=config.geometry.span_m;U=Re*config.flow.mu_Pa_s/(config.flow.rho_kg_m3*chord);R=[cos(alpha),sin(alpha);-sin(alpha),cos(alpha)];
[eta,w]=fem3d.gauss01(8);nq=size(rows,1)*numel(eta);I=cell(nq,1);J=I;V=I;points=zeros(nq,3);loads=points;count=0;
for j=1:numel(eta),for k=1:size(rows,1)
 count=count+1;xy=[rows(k,2)*chord,(eta(j)-.5)*span];hit=cad.ray(mesh,xy);assert(numel(hit.z)>=2,'cad:Traction','Missing material skin.');if rows(k,1)==1,face=hit.face(end);else,face=hit.face(1);end
 ids=mesh.faces(face,:);[point,N,~,ok]=cad.face_point(mesh,ids,xy);assert(ok);I{count}=repmat(count,numel(ids),1);J{count}=ids';V{count}=N';points(count,:)=point;
 force=rows(k,8:9)*R;loads(count,:)=config.flow.rho_kg_m3*U^2*chord*span*w(j)*[force(1),0,force(2)];
end,end
H=sparse(vertcat(I{:}),vertcat(J{:}),vertcat(V{:}),nq,size(mesh.nodes,1));F=H'*loads;
transfer=struct('H',H,'points',points,'pointForces',loads,'nodalForces',F,'speed_ms',U,'Re',Re,'forceDefect',norm(sum(F,1)-sum(loads,1)),'momentDefect',norm(sum(cross(mesh.nodes,F,2),1)-sum(cross(points,loads,2),1)),'scope','span-uniform extrusion of computed low-Re 2D viscous skin traction; original upper/lower material FE interpolation');
end
