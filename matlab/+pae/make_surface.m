function surface = make_surface(planform, nx, ny)
%MAKE_SURFACE Mean aerodynamic surface. All lengths in metres, angles rad.
% planform.spanStations: increasing CAD-independent span coordinates.
% planform.leadingEdge: N-by-3; its second coordinate equals spanStations.
% planform.chord: N-by-1. Optional twist (positive nose-up), camberFcn(xi,s).
% Surface represents the ENTIRE finite lifting surface, not a symmetry half.
arguments
    planform struct
    nx (1,1) double {mustBeInteger,mustBePositive} = 6
    ny (1,1) double {mustBeInteger,mustBePositive} = 16
end
s = planform.spanStations(:);
le = planform.leadingEdge;
c = planform.chord(:);
assert(numel(s)>=2 && all(diff(s)>0),'pae:Planform','Span stations must increase.');
assert(isequal(size(le),[numel(s),3]) && numel(c)==numel(s), ...
    'pae:Planform','Invalid leadingEdge/chord dimensions.');
assert(all(isfinite([s;le(:);c])) && all(c>0),'pae:Planform','Invalid planform values.');
spanRange = s(end)-s(1);
assert(max(abs(le(:,2)-s))<1e-10*max(1,spanRange), ...
    'pae:Planform','leadingEdge(:,2) must equal spanStations.');
xi = (1-cos(linspace(0,pi,nx+1)))/2;
eta = (1-cos(linspace(0,pi,ny+1)))/2;
span = s(1)+spanRange*eta;
LE = interp1(s,le,span,'linear');
chord = interp1(s,c,span,'linear');
twist = zeros(size(span));
if isfield(planform,'twist')
    assert(numel(planform.twist)==numel(s) && all(isfinite(planform.twist)), ...
        'pae:Planform','Invalid twist values.');
    twist = interp1(s,planform.twist(:),span,'linear');
end
[XI,ETA] = ndgrid(xi,eta);
camber = zeros(size(XI));
if isfield(planform,'camberFcn')
    camber = planform.camberFcn(XI,ETA); % dimensionless camber z/c
    assert(isequal(size(camber),size(XI)) && all(isfinite(camber(:))), ...
        'pae:Planform','camberFcn must return a finite grid of z/c.');
end
nodes = zeros(nx+1,ny+1,3);
for j=1:ny+1
    xlocal = chord(j)*xi(:);
    zlocal = chord(j)*camber(:,j);
    nodes(:,j,1) = LE(j,1)+xlocal*cos(twist(j))+zlocal*sin(twist(j));
    nodes(:,j,2) = LE(j,2);
    nodes(:,j,3) = LE(j,3)-xlocal*sin(twist(j))+zlocal*cos(twist(j));
end
surface = struct('nodes',nodes,'xi',XI,'eta',ETA,'nx',nx,'ny',ny, ...
    'referenceArea',trapz(s,c),'referenceChord',trapz(s,c.^2)/trapz(s,c), ...
    'span',spanRange,'source','parametric mean surface');
end
