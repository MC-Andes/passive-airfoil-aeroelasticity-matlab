function surface = surface(mesh,config,nx,ny)
%SURFACE Material outer-skin mean surface with explicit upper/lower FE operators.
% Positive/negative pressures are not reconstructed from a thin-surface load.
% The first/last 0.1% chord is excluded to avoid singular edge inversion.
xi=linspace(.001,.999,nx+1); eta=(1-cos(linspace(0,pi,ny+1)))/2;
nn=size(mesh.nodes,1); np=(nx+1)*(ny+1); I={}; J={}; V={}; upper=zeros(np,3); lower=upper; count=0;
for j=1:ny+1
    y=(eta(j)-.5)*config.geometry.span_m;
    for i=1:nx+1
        count=count+1; x=xi(i)*config.geometry.chord_m; hit=cad.ray(mesh,[x,y]);
        assert(numel(hit.z)>=2,'cad:Surface','Outer envelope missing at x/c=%g, eta=%g.',xi(i),eta(j));
        for side=1:2
            if side==1,f=hit.face(end);else,f=hit.face(1);end
            face=mesh.faces(f,:); xe=mesh.nodes(face,:);
            [~,w,~,ok]=cad.face_point(mesh,face,[x,y]);
            assert(ok,'cad:Surface','Cannot invert outer material skin.');
            I{end+1}=repmat(count,numel(face),1); J{end+1}=face'; V{end+1}=w'/2; %#ok<AGROW>
            if side==1,upper(count,:)=w*xe;else,lower(count,:)=w*xe;end
        end
    end
end
H=sparse(vertcat(I{:}),vertcat(J{:}),vertcat(V{:}),np,nn); points=H*mesh.nodes;
assert(norm(points-(upper+lower)/2,inf)<1e-10,'cad:Surface','Mean-skin interpolation mismatch.');
[XI,ETA]=ndgrid(xi,eta);
surface=struct('nodes',reshape(points,nx+1,ny+1,3),'H',H,'Oxyz',kron(H,speye(3)), ...
    'upperPoints',upper,'lowerPoints',lower,'nx',nx,'ny',ny,'xi',XI,'eta',ETA, ...
    'referenceArea',config.geometry.chord_m*config.geometry.span_m, ...
    'referenceChord',config.geometry.chord_m,'span',config.geometry.span_m, ...
    'interpretation','CAD outer-skin mean; displacement is the average of exposed upper/lower material points', ...
    'loadPartition','work-conjugate thin-surface forces split equally to outer skins; no resolved Cp or cavity pressure');
end
