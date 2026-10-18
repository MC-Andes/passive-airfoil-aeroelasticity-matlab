function sensor = observers(mesh,c)
%OBSERVERS Material TE points in positive span coordinates of the updated CAD.
eta=c.sensors.eta(:);x=min(mesh.nodes(:,1))+c.sensors.x_fraction_of_actual_cad_extent*range(mesh.nodes(:,1));
I=[];J=[];V=[];points=zeros(numel(eta),3);
for k=1:numel(eta)
 hit=cad.ray(mesh,[x,eta(k)*c.geometry.span_m]);assert(~isempty(hit.z),'morphing:Sensor','No material skin intersection.');
 face=mesh.faces(hit.face(end),:);[p,N,~,ok]=cad.face_point(mesh,face,[x,eta(k)*c.geometry.span_m]);assert(ok);
 I=[I;repmat(k,numel(face),1)];J=[J;face'];V=[V;N'];points(k,:)=p; %#ok<AGROW>
end
H=sparse(I,J,V,numel(eta),size(mesh.nodes,1));sensor=struct('H',H,'O',kron(H,sparse([0,0,1])),'points',points,'eta',eta,'x_over_nominal_c',x/c.geometry.chord_m,'scope','actual finite-CAD TE offset, positive root-to-tip span, upper-skin material sensors');
end
