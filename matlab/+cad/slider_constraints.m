function mechanism = slider_constraints(mesh,config,scenario)
%SLIDER_CONSTRAINTS Explicit hypothesis from the cross-section's internal bar.
% CAD faces 7 and 9 bound parallel horizontal bar faces around the spanwise
% hole. They define a CANDIDATE slider body, not proof of actual guide contact.
% A: connected elastic solid. B: block its six mean rigid coordinates.
% C: guide its mean motion along chord x; elasticity/connectivity preserved.
labels=cad.face_labels(mesh);selected=mesh.faces(ismember(labels,[7,9]),:);region=unique(selected);
assert(numel(region)>=6,'cad:Slider','Candidate internal-bar surfaces are absent.');
quadrature=fem3d.surface_quadrature(mesh,selected,3);X=quadrature.points;
center=sum(quadrature.weights.*X,1)/sum(quadrature.weights);local=X-center;R=zeros(3*size(X,1),6);
R(1:3:end,1)=1;R(2:3:end,2)=1;R(3:3:end,3)=1;
for k=1:3,w=zeros(1,3);w(k)=1;U=cross(repmat(w,size(X,1),1),local,2);R(:,k+3)=reshape(U',[],1);end
weighted=repelem(quadrature.weights,3).*R;
observation=sparse((R'*weighted)\(weighted'*kron(quadrature.H,speye(3))));
if strcmp(scenario,'A'),C=sparse(0,size(observation,2));
elseif strcmp(scenario,'B'),C=observation;
elseif strcmp(scenario,'C'),C=observation(2:6,:);
else,error('cad:Slider','Unknown slider scenario.');end
mechanism=struct('C',C,'d',zeros(size(C,1),1),'rigidObservation',observation, ...
    'coordinateDefinition','area-weighted least-squares rigid fit on physical CAD faces 7/9','nodes',region,'center',center,'scenario',scenario,'guideDirection',[1,0,0], ...
    'provenance','assumed mean-coordinate guide from parallel CAD faces 7/9; authorized geometric/literature hypothesis', ...
    'limits','no nodes separated; not identification of actual joints, preload or friction');
end
