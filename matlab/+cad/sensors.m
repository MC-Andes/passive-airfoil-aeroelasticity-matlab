function sensor = sensors(mesh,config)
%SENSORS Fixed material near-TE stations on the exposed upper skin.
% Same canonical x/y on every mesh; no nearest-node snapping. H is immutable
% during deformation. The finite x/c offset is reported, not called the sharp edge.
eta=config.sensors.eta(:); count=numel(eta); nn=size(mesh.nodes,1);
x=config.geometry.chord_m*config.sensors.x_over_c;
y=(eta-.5)*config.geometry.span_m; I=[]; J=[]; V=[]; points=zeros(count,3); faces=zeros(count,1);
for k=1:count
    h=cad.ray(mesh,[x,y(k)]);
    assert(~isempty(h.z),'cad:Sensor','Near-TE material point has no surface intersection.');
    f=h.face(end); face=mesh.faces(f,:); xe=mesh.nodes(face,:);
    [~,weights,~,ok]=cad.face_point(mesh,face,[x,y(k)]);
    assert(ok,'cad:Sensor','Cannot invert material sensor face.');
    I=[I;repmat(k,numel(face),1)]; J=[J;face']; V=[V;weights']; %#ok<AGROW>
    points(k,:)=weights*xe; faces(k)=f;
end
H=sparse(I,J,V,count,nn); O=kron(H,sparse([0,0,1]));
sensor=struct('H',H,'O',O,'points',points,'eta',eta,'faces',faces, ...
    'x_over_c',config.sensors.x_over_c,'direction',[0,0,1], ...
    'interpretation','signed z displacement of fixed material points near TE; finite chord offset');
end
