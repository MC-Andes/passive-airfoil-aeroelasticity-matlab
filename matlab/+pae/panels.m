function panel = panels(nodes)
%PANELS Geometry and interpolation operators for a horseshoe VLM.
[nxc,nyc,ncomp] = size(nodes);
assert(nxc>=2 && nyc>=2 && ncomp==3 && all(isfinite(nodes(:))), ...
    'pae:Surface','Nodes must be a finite (nx+1)-by-(ny+1)-by-3 grid.');
points = reshape(nodes,[],3);
ids = reshape(1:nxc*nyc,nxc,nyc);
corners = [reshape(ids(1:end-1,1:end-1),[],1), ...
    reshape(ids(2:end,1:end-1),[],1), ...
    reshape(ids(1:end-1,2:end),[],1), ...
    reshape(ids(2:end,2:end),[],1)];
np = size(corners,1); nn = size(points,1);
rows = repmat((1:np)',1,4);
panel.Hforce = sparse(rows(:),corners(:),repmat([.375,.125,.375,.125],np,1),np,nn);
panel.Hcontrol = sparse(rows(:),corners(:),repmat([.125,.375,.125,.375],np,1),np,nn);
p00 = points(corners(:,1),:); p10 = points(corners(:,2),:);
p01 = points(corners(:,3),:); p11 = points(corners(:,4),:);
panel.A = .75*p00+.25*p10;
panel.B = .75*p01+.25*p11;
panel.control = panel.Hcontrol*points;
panel.forcePoint = panel.Hforce*points;
chord = .5*(p10-p00+p11-p01);
span = .5*(p01-p00+p11-p10);
normal = cross(chord,span,2);
area = vecnorm(normal,2,2);
assert(all(area>1e-14),'pae:Surface','Collapsed panels detected.');
panel.normal = normal./area;
panel.area = area;
panel.corners = corners;
end
